import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import matplotlib.font_manager as fm
from matplotlib.legend_handler import HandlerPatch
from matplotlib.path import Path
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import xarray as xr
import pandas as pd
import seaborn as sns
import seaborn_image as isns
import cmasher as cm
import warnings
import string
from scipy.stats import linregress
import numpy as np
import copy
import glob
import os
import math
import xesmf as xe


warnings.filterwarnings('ignore')


fm.FontProperties(family='serif', size=14)
isns.set_context(mode="notebook", fontfamily="serif")


def create_discrete_cmap(cmap, n):

    colors = cmap(np.linspace(0, 1, n))
    discrete_cmap = mcolors.ListedColormap(colors)
    bounds = np.linspace(0, 100, n+1)
    norm = mcolors.BoundaryNorm(bounds, discrete_cmap.N)
    
    return discrete_cmap, norm
    
n = 25


def month_x_labels(ax):
    
    month_label = ['J','F','M','A','M','J','J','A','S','O','N','D']
    ax.set_xticks(np.arange(0, 12, 1))
    ax.set_xticklabels(month_label, fontsize=12)


def var_units():
    # variable dictionary with units
    var_units = {
                        'ALK' : 	{'Alkalinity' : r'μ mol $\mathregular{kg^{-1}}$'},
                        'DIC' : 	{'DIC' : r'μ mol $\mathregular{kg^{-1}}$'},
                        'co2flux' : 	{'$\mathregular{CO_{2}}$ flux' : r'kg $\mathregular{m^{-2} \ y^{-1}}$'},
                        'fco2' : 	{'Ocean p$\mathregular{CO_{2}}$' : r'μatm'},
                        'ph': 		{'pH': 'unit'},
                        'somxl010': 	{'MLD': 'm'},
                        'bgc_diag_pp':	{'NPP': r'mmol C $\mathregular{m^{-3} \ y^{-1}}$'},
                        'sosstsst': 	{'SST': r'$^\circ$C'},                   
                }

    return var_units


def lat_lon_labels():
    
    # degree symbol in string
    d = '$^\circ$' 
    # latitude coordinates
    lat_labels = [f'38{d}N', f'44{d}N',f'50{d}N', f'56{d}N', f'62{d}N', f'69{d}N',f'76{d}N'] 
    # longitude coordinates
    lon_labels = [f'28{d}W', f'20{d}W', f'12{d}W',f'4{d}W', f'4{d}E',f'12{d}E']

    return lat_labels, lon_labels


def seasonal_avg(data, var):

    # days in a month
    m_length = data.time_counter.dt.days_in_month
    
    # weighted monthly average
    monthly_avg = (data[var] * m_length).groupby("time_counter.month").mean(dim="time_counter") / m_length.groupby("time_counter.month").mean(dim="time_counter")
    monthly_avg = monthly_avg.to_dataset(name=var)
    
    return monthly_avg


def ocean_mask(data_path, var): 
    # ocean mask

    mesh_mask = xr.open_dataset(data_path+'mesh_mask.nc')
    dx, dy = mesh_mask.e1t, mesh_mask.e2t
    grid_cell_area = (dx * dy).isel(t=0)

    mesh_mask['tmask'] = mesh_mask['tmask'].where(mesh_mask['tmask']!=0)
    mesh_mask['tmask'] = mesh_mask['tmask'].isel(t=0,z=0)*grid_cell_area
    mesh_mask['tmask'] = mesh_mask['tmask'].where(mesh_mask['tmask'])

    if var in ['sosstsst','bgc_diag_pp']:
        europe_mask = mesh_mask.sel(x=slice(519, 595), y=slice(329, 435)) # m-2
    else:
        europe_mask = mesh_mask.sel(x=slice(520, 595), y=slice(330, 435)) # m-2
    
    return europe_mask


def coastline_mask(data_path):
    # alkalkinity mask

    alk_mask = xr.open_dataset(data_path+'alkalinity_mask_y2100.nc')
    alk_mask = alk_mask.where(alk_mask)
    alk_mask = alk_mask.sel(x=slice(520, 595), y=slice(330, 435))
    alk_mask = xr.where(alk_mask.notnull(), 1, np.nan)
    alk_mask = alk_mask.where(alk_mask)
    
    return alk_mask


degree = '$^\circ$'
data_path = '/Users/chiaraciscato/Desktop/GEOMAR/2025_ciscato_oae_seasonality/datasets/last_decade/'


def fineline():
    
    for location in ['left','bottom', 'right', 'top']:
        ax.spines[location].set_linewidth(0.4)
