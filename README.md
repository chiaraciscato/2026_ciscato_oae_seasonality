## Impacts of Simulated Coastal Ocean Alkalinity Enhancement on the Seasonal Carbon Cycle in European Waters under a Low- and a High-Emission Scenario

This repository accompanies this paper : Ciscato, C., Mehendale, N., Kemena, T., Avrutin, S., Keller, D. P. 2026 (Earth System Dynamics)

### Overview 

A seasonal analysis was performed on different carbon cycle variables, with the objective of defining the impacts of Ocean Alkalinity Enhancement on their monthly (seasonal) cycle. Two reference scenarios are used: SSP1-2.6 (low warming) and SSP3-7.0 (high warming). The model domain is the European coastline (excluding the Mediterranean and the Baltic seas).

### Repository content

This repository contains three folders:

- the [scripts](scripts) folder contains the python scripts that process the datasets and plot the final figures. 
- the [out](out) folder contains the final figures conceived to showcase the results.

### To keep in mind

The python scripts assume that the input data have already been sliced to the region of focus, namely the European coastline. The exact extremes are noted in the script that addresses the seasonal analysis of the variables. Additionally, as values over land are set to zero, a mask should be applied to all variables to exclude land values from the calculations. 

## Author

- [@chiaraciscato](https://github.com/chiaraciscato)

