script_file = "general_functions.py"
with open(script_file) as f:
    content = f.read()

exec(content)

# variable names
vars = ['ALK', 'DIC', 'fco2', 'co2flux', 'ph', 'mld', 'npp', 'sst']
# scenario names
scenarios = ['baseline_126', 'oae_126', 'baseline_370', 'oae_370']

file_path = {}

# open datasets
for var in vars:
    file_path[var] = {}

    for s in scenarios:
        file_path[var][s] = data_path + var + f'_{s}.nc'
        
# variable units
units = var_units()

# latitude, longitude coordinates
lat_labels, lon_labels = lat_lon_labels()


def vertical_weights(data, var, depth_dataset, oper):
    
    depth = data['deptht']

    depth_values = depth.values
    layer_thickness = np.empty_like(depth_values)
    layer_thickness[:-1] = np.diff(depth_values)
    layer_thickness[-1] = layer_thickness[-2] 

    layer_thickness = xr.DataArray(
        layer_thickness,
        dims=['deptht'],
        coords={'deptht': depth}
    )

    if oper == 'ave':
        weighted_sum = (data[var] * layer_thickness).sum(dim=depth, skipna=True)
        result = (weighted_sum / depth_dataset).to_dataset(name=var)

    elif oper == 'sum':
        result = ((data[var] * layer_thickness).sum(dim=depth, skipna=True)).to_dataset(name=var)

    return result


def horizontal_weights(data, var, mask, mask_var, oper):

    if oper == 'ave':
        weighted_sum = (data[var] * mask[mask_var]).sum(dim=['x', 'y'], skipna=True)
        total_weight = mask[mask_var].sum(dim=['x', 'y'], skipna=True)
        result = weighted_sum / total_weight

    elif oper == 'sum':
        result = (data[var] * mask[mask_var]).sum(dim=['x', 'y'], skipna=True)

    result = result.to_dataset(name=var)
    
    return result


def fig_plot_map(data, var, ax, vmin, vmax, oae_index, base_index):
    # get the difference between oae and baseline over time, then compute the 10-year mean
    if var != 'DIC':
        da = ((list(data.values())[oae_index]-list(data.values())[base_index]).mean('year'))[var]
    else:
        da = (list(data.values())[oae_index].isel(year=-1)-list(data.values())[base_index].isel(year=-1))[var]
    
    im = ax.pcolormesh(
        da,
        cmap=cmap,
        vmin=vmin,vmax=vmax,
        shading='auto')

    return im


def process_scenario(var, file_path):
    # load datasets
    var_datasets = {
        scenario: xr.open_dataset(path, decode_times=True)
        for scenario, path in file_path[var].items()
    }
    # load mixed layer depth datasets
    mld_datasets = {
        scenario: xr.open_dataset(path, decode_times=True)
        for scenario, path in file_path['mld'].items()
    }
    # start dictionaries
    nodepth_datasets = {}
    map_datasets = {}
    integral = {}
    
    for key, d in var_datasets.items():
        if 'deptht' in d.dims:
            depth = 'deptht'
            max_depth = d[depth].isel({depth:d[var].argmax(dim=depth)})
        
        d = d.where(d).sel(time_counter=slice('2090','2100'))

        d_mld = mld_datasets[key]

        # if var in ['ALK','DIC']:
        #     # mask = (data.deptht > 200) 
            
        #     mask = (dataset['deptht'] < d_mld['somxl010'])
        #     dataset[var] = dataset[var].where(mask)
        #     dataset = vertical_weights(dataset, var, d_mld, 'somxl010', 'ave')
            
        if var in ['ALK']:
            d = d.isel({depth:0})
        
        elif var in ['DIC']: 
            ocn_mask = ocean_mask(data_path, var)

            intgr = vertical_weights(d, var, max_depth, 'sum')
            # intgr = horizontal_weights(intgr, var, ocn_mask, 'tmask', 'ave')
            
            intgr[var] = intgr[var] * 1025 * 1e-6 * 12.01 * 1e-3 # C kg m-3         
            intgr = intgr.groupby('time_counter.year').mean('time_counter')
            integral[key] = intgr        
            
        elif var in ['co2flux']:
                                          
            d[var] = d[var] * 31536000
    
        nodepth_datasets[key] = d

        d = (d.groupby('time_counter.year').max(dim='time_counter'))-(d.groupby('time_counter.year').min(dim='time_counter'))
        map_datasets[key] = d

    return nodepth_datasets, map_datasets, integral

# space weighting 
def space_weighting(processed_datasets, var):

    european = {}
    coastline = {}

    for key, d in processed_datasets.items():
        
        # european average
        ocn_mask = ocean_mask(data_path, var)
        eu = horizontal_weights(d, var, ocn_mask, 'tmask', 'ave')
        european[key] = eu

        # coastline average
        cst_mask = coastline_mask(data_path + '../alk_mask/')
        coast = horizontal_weights(d, var, cst_mask, 'alk_flux', 'ave')
        coastline[key] = coast

    return european, coastline

# calculate weighted average and seasonal amplitudes
def seasonal_amplitude(scenarios, var):
    year_scenario = {key: seasonal_avg(d, var) for key, d in scenarios.items()}
    
    delta_scenario = {
        f'delta_{baseline_key}': year_scenario[baseline_key] - year_scenario[oae_key]
        for baseline_key in year_scenario
        # zip base and oae scenarios
        if baseline_key.startswith('baseline') and (oae_key := baseline_key.replace('baseline', 'oae')) in year_scenario
    }

    return year_scenario, delta_scenario

def data_processing(var, file_path):

    print('processing  ' + var)

    nodepth_datasets, map_amplitude, integral = process_scenario(var, file_path)
    
    european, coastline = space_weighting(nodepth_datasets, var)
    final_eu = seasonal_amplitude(european, var)
    final_point = seasonal_amplitude(coastline, var)

    if var in ['co2flux','DIC']:
        return final_eu, final_point, map_amplitude, integral
    else:
        return final_eu, final_point, map_amplitude

    print('finished ' + var)




eu_alk, ps_alk, ampl_alk = data_processing('ALK', file_path)
eu_pco2, ps_pco2, ampl_pco2 = data_processing('fco2', file_path)
eu_co2flux, ps_co2flux, ampl_co2flux, co2flux_integral = data_processing('co2flux', file_path)

# eu_dic, ps_dic, ampl_dic, dic_integral = data_processing('DIC', file_path)
# eu_ph, ps_ph, ampl_ph = data_processing('ph', file_path)
# eu_mld, ps_mld, ampl_mld = data_processing('mld', file_path)


print('generating lineplot')

s, m, l = 16, 20, 23

datasets = [eu_alk, eu_pco2, eu_co2flux, ps_alk, ps_pco2, ps_co2flux]
vars = ['ALK','fco2','co2flux','ALK','fco2', 'co2flux']

# datasets = [eu_dic, eu_ph, eu_mld, ps_dic, ps_ph, ps_mld]
# vars = ['DIC','ph','somxl010','DIC','ph', 'somxl010']

labels = ['Baseline', 'OAE']
colors = ['lightcoral', 'mediumturquoise']

def lineplot_vars(ds, var, ax1, ax2, ax3, ax4, labels, colors):
    # abs 
    ax1.plot(list(ds[0].values())[0][var], label=labels[0], color=colors[0], linewidth=4) # baseline ssp126
    ax1.plot(list(ds[0].values())[1][var], label=labels[1], color=colors[1], linewidth=4) # oae ssp126
    ax3.plot(list(ds[0].values())[2][var], label=labels[0], color=colors[0], linewidth=4) # baseline ssp370
    ax3.plot(list(ds[0].values())[3][var], label=labels[1], color=colors[1], linewidth=4) # oae ssp370

    # delta
    ax2.plot(list(ds[1].values())[0][var], color='darkseagreen', linewidth=2) # delta ssp126
    ax4.plot(list(ds[1].values())[1][var], color='darkseagreen', linewidth=2) # delta ssp370

    ax1.set_ylabel(f"[ {list(units[var].values())[0]} ]", fontsize=s, labelpad=s)
    ax2.set_ylabel(f"Δ", fontsize=s, labelpad=s)

    ax1.set_title('SSP1-2.6', fontsize=m, y=1.02)
    ax3.set_title("SSP3-7.0", fontsize=m, y=1.02)

    ax3.yaxis.set_tick_params(labelleft=False)
    ax4.yaxis.set_tick_params(labelleft=False)

    if var in ['co2flux']:
        ax1.axhline(0, color='black', linewidth=0.5)
        ax3.axhline(0, color='black', linewidth=0.5)

gs0 = gridspec.GridSpec(6, 3, hspace=0.2)

axes = []

fig = plt.figure(figsize=(30, 48)) 

for i in range(6):
    gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[i], height_ratios=[1, 0.4], width_ratios=[1, 1], hspace=0, wspace=0)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax4 = fig.add_subplot(gs[1, 1], sharey=ax2)

    axes.extend([ax1, ax2, ax3, ax4])
    
    lineplot_vars(datasets[i], vars[i], ax1, ax2, ax3, ax4, labels, colors)

    if i == 0:
            ax1.legend(labels=labels, loc='best', frameon=False, handlelength=.4)

for i, ax in enumerate(axes):
    month_x_labels(ax)
    ax.tick_params(axis='both', which='major', labelsize=s)
    fineline()

for n, ax in enumerate([axes[i] for i in range(0, 24, 2)]):
    ax.text(0.87, 0.4, f"{string.punctuation[7]}{string.ascii_lowercase[n]}{string.punctuation[8]}",
            transform=ax.transAxes, size=m)

column_titles = [f"Surface alkalinity", f"{list(units[vars[1]].keys())[0]}", f"{list(units[vars[2]].keys())[0]}"]
for i, title in enumerate(column_titles):
    fig.text(0.195 + i * 0.283, 0.895, title, fontsize=l)

row_titles = ['European Avg.', 'Coastline Avg.']
for i, title in enumerate(row_titles):
    fig.text(0.06, 0.8 - i * 0.13,title, ha='left', fontsize=l, rotation=90)

# plt.savefig('/Users/chiaraciscato/Desktop/GEOMAR/2025_ciscato_oae_seasonality/out/figure_3.png')
# plt.savefig('/Users/chiaraciscato/Desktop/GEOMAR/2025_ciscato_oae_seasonality/out/lineplot_supplement.png', transparent=True)

pass


print('generating maps')

cmap, norm = create_discrete_cmap(cm.seasons, 40)
cmap.set_bad(color='silver')


fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, figsize=(30, 22))
plt.subplots_adjust(wspace=0.11, hspace=0.09)

plt.rcParams['axes.linewidth'] = .2

datasets = [ampl_alk, ampl_pco2, ampl_co2flux]
vars = ['ALK','fco2','co2flux','ALK','fco2','co2flux']
num1,num2,num3 = 400,200,1

# datasets = [ampl_dic, ampl_ph, ampl_mld]
# vars = ['DIC','ph','somxl010','DIC','ph', 'somxl010']
# num1,num2,num3 = 600,.1,200

axes = [ax1, ax2, ax3, ax4, ax5, ax6]

im1 = fig_plot_map(datasets[0], vars[0], ax1, vmin=-num1, vmax=num1, oae_index=1, base_index=0)
im2 = fig_plot_map(datasets[1], vars[1], ax2, vmin=-num2, vmax=num2, oae_index=1, base_index=0)
im3 = fig_plot_map(datasets[2], vars[2], ax3, vmin=-num3, vmax=num3, oae_index=1, base_index=0)
im4 = fig_plot_map(datasets[0], vars[0], ax4, vmin=-num1, vmax=num1, oae_index=3, base_index=2)
im5 = fig_plot_map(datasets[1], vars[1], ax5, vmin=-num2, vmax=num2, oae_index=3, base_index=2)
im6 = fig_plot_map(datasets[2], vars[2], ax6, vmin=-num3, vmax=num3, oae_index=3, base_index=2)

scenarios = ['SSP1-2.6','SSP3-7.0']

for i, (ax, scen) in enumerate(zip([ax1, ax4], scenarios)):
    ax.set_ylabel(f'{scenarios[i]}', labelpad=14)
    
for ax in axes: 
    fineline()
    lat_labels = lat_labels
    lon_labels = lon_labels
        
    if ax in [ax1,ax4]:
        ax.set_yticklabels(lat_labels, fontsize=s)        
    else:
        ax.yaxis.set_tick_params(labelleft=False)
        
    if ax in [ax4,ax5,ax6]:
        ax.set_xticks([0,15,30,45,60,75])
        ax.set_xticklabels(lon_labels, fontsize=s)
    else:
        ax.xaxis.set_tick_params(labelleft=False) 
        
    ax1.set_title(f"MLD-avg. {list(units[vars[0]].keys())[0]} [ {list(units[vars[0]].values())[0]} ]", fontsize=m, y=1.22)
    ax2.set_title(f"Surface {list(units[vars[1]].keys())[0]} [ {list(units[vars[1]].values())[0]} ]", fontsize=m, y=1.22),
    ax3.set_title(f"{list(units[vars[2]].keys())[0]} [ {list(units[vars[2]].values())[0]} ]", fontsize=m, y=1.22),
        
for n, ax in enumerate(axes):
    ax.text(0.88, 0.9, string.punctuation[7] + string.ascii_lowercase[n] + string.punctuation[8],
            transform=ax.transAxes, size=m, color='black')

ims = [im1,im2,im3]
axes = [[ax1,ax4],[ax2,ax5],[ax3,ax6]]

for im, axes in zip(ims, axes):
    cbar = fig.colorbar(im, ax=axes, aspect=15, pad=0.02, location='top', extend='both')

pass

# plt.savefig('/Users/chiaraciscato/Desktop/GEOMAR/2025_ciscato_oae_seasonality/out/map_updated.png')
# plt.savefig('/Users/chiaraciscato/Desktop/GEOMAR/2025_ciscato_oae_seasonality/out/map_supplement.png', transparent=True)
