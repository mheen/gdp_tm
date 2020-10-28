# gdp_tm
Deriving transport matrices from observed drifter locations from NOAA's Global Drifter Program

## Getting started
Head to http://www.aoml.noaa.gov/envids/gld/FtpInterpolatedInstructions.php and download NOAA's Global Drifter Program (GDP) drifter data.
You'll need both the data (buoydata_xxx.dat) and the meta data (dirfl_xxx.dat) files.

In the gdp_tm/drifters folder create a file "dirs.json". This file points to relevant file locations that "read_gdp_data" uses.
It should contain the following fields:
```
{
 "dirs":
  {
   "input_dir": "my_directory/with_gdp/original_data/",
   "output_path": "my_directory/with_gdp/matlab_data.mat",
   "data_filenames": ["buoydata_1_5000.dat","buoydata_5001_10000.dat"],
   "meta_data_filenames": ["dirfl_1_5000.dat","dirfl_5001_10000.dat"]
  }
 }
 ```
 
 Now run "gdp_data = read_gdp_data()". This will load all the GDP drifter data into a Matlab data structure
 (and save it to the .mat file you specified in "dirs.json" so you don't have to do this again all the time).
 
## Creating transport matrices and running transport matrix simulations
To run a transport matrix simulation, you need to:
1. Load the .mat file with the GDP data created by "read_gdp_data".
2. Create a .ini file containing your model configuration settings (see the section "Model configuration" below).

After that, run a simulation by:
```
config = ModelConfiguration('config.ini');
tm = TransportMatrix(gdp_data,config);
simulation = TmModel(config,tm);
```
Check out the script "tm/main.m" that does this, and plots some of the simulation results as well.

The sections below explain these steps in more detail.

### Model configuration
To derive transport matrices and run simulations, you first need to create a model configuration using "ModelConfiguration".
You can create a model configuration in two ways:
1. Directly pass settings to ModelConfiguration, for example:
   ```
   config = ModelConfiguration('drogued_status','undrogued');
   ```
   To load default settings, don't pass any arguments to ModelConfiguration.
2. Pass a config.ini file to ModelConfiguration, for example:
   ```
   config = ModelConfiguration('config.ini');
   ```
   with for example the config.ini file:
   ```
   [tm model configuration]
   drogued_status=undrogued
   ```
   
All default values and an explanation of all the settings are given below.
```
--------------------------------------------------------------------------
Default settings
--------------------------------------------------------------------------
drogued_status=all
dx=1
dt=60
run_time=50
output_interval=1
sinks=0
normalise=0
initial.type=global_uniform
initial.lon=[]
initial.lat=[]
initial.tracer_per_grid_cell=1
initial.years_to_add_sources=[]
days_in_year=360
specify_time.type=all
specify_time.times=[]
   
--------------------------------------------------------------------------
Additional options/information
--------------------------------------------------------------------------
- drogued_status [string]: "all", "drogued" or "undrogued"

- dx [degrees]: Grid cell size of the transport matrix.

- dt [days]: Separation time between drifter pairs in transport matrix
   construction.

- run_time [years]: Simulation run time.

- output_interval [integer]: Interval for which output is saved.
   1 indicates every time step, 2 every other time step, etc.

- sinks [logical]: Indicates whether sinks should be fixed in the transport
   matrix (1) or removed (0).

- normalise [logical]: Indicates whether simulated tracer output should
  be normalised or not (1 uses normalisation, 0 does not).
  Normalisation is done by multiplying the amount of tracer in a specific
  grid cell c(x,y) by the number of total grid cells N divided by the
  total amount of tracer C in the simulation at that time,
  TAF = c(x,y)*N/C.
  TAF is the tracer amplification factor as defined by van Sebille
  et al. (2012): "If, for example, in some grid box TAF(x,y)=20, then
  twenty times more tracer is found within that grid box than if all
  tracer is uniformly distributed over the global ocean."
  Normalisation is always done if the initial condition injects plastic
  sources, such as when using initial.type=jambeck2015. In this case, if
  normalisation is set to 0 here, this is overruled when the
  ModelConfiguration object is constructed.

- initial: Consists of several subfields:
  * type [string]: "global_uniform": Uniformly distributed tracer over
                    global oceans.
                   "uniform": Uniformly distributed tracer in specified
                    longitude and latitude range.
                    "jambeck2015": Releases tracer from coasts according
                    to Jambeck et al. (2015) distribution.
                    Note: this functionality has not yet been built.
                    "lebreton2017": Releases tracer from rivers according
                    to Lebreton et al. (2017) distribution.
                    Note: this functionality has not yet been built.

  To use with type "global_uniform" and "uniform":
  * tracer_per_grid_cell [real number]: Amount of tracer initially
    released in a grid cell.

  To use with type "uniform":
  * lon [array]: Longitude range, for example: [19,117].
  * lat [array]: Latitutde range, for example: [-45,-15].

  To use with type "jambeck2015" and "lebreton2017":
  * years_to_add_source [years]: Specifies how long source term
    (initial condition) should be added to simulation.

- days_in_year [integer]: Must be a multiple of the separation time dt.
  This is used to determine the number of model iterations needed to run
  a simulation for the duration of time specified by run_time.
  In the default settings, dt=60 days and days_in_year is therefore
  rounded down from 365 days to 360 days. In this case, 6 model
  iterations simulate 1 year.

- specify_time: Consists of several subfields:
  * type [string]: "all": Use all available drifter times to construct
                   the transport matrix.
                   "month": Use only drifter times from specific months
                   to construct the transport matrix.
                   "year": Use only drifter times from specific years to
                   construct the transport matrix.

  To use with type "month" or "year":
  * times [array]: Array with relevant months or years that will be used
                   to construct the transport matrix.
                   
```

### Transport matrix
Create a transport matrix by:
```
tm = TransportMatrix(gdp_data,config);
```
"TransportMatrix" derives a (non-normalised) transport matrix from the GDP drifter data, using the settings in a "ModelConfiguration". All sink locations are then iteratively removed from the transport matrix. Finally, the transport matrix is normalised. The normalised transport matrix ("tmn") is used for simulations.

Information about the derivation of the transport matrix is also stored. This includes:
* locations of sinks,
* locations of undefined grid cells (no drifter measurements in these cells),
* the number of individual drifters used to construct the transport matrix,
* the number of observed drifter locations used to construct the transport matrix,
* the number of observed drifter locations in each month,
* the number of observed drifter locations in each year.

The final transport matrix is saved to a .mat file in an "output/" folder. If a transport matrix with the requested model configuration already exists when calling TransportMatrix, it is loaded from this .mat file instead of deriving it again.

### Transport matrix simulations
Run a transport matrix simulation by:
```
simulation = TmModel(config,tm);
```
"TmModel" creates an initial condition and then runs a simulation, using the settings in a "ModelConfiguration".

The simulation result is saved to a .mat file in an "output/" folder. If a simulation with the requested model configuration already exists when calling TmModel, it is loaded from this .mat file instead of run again.

## Postprocessing and plotting
Under "tm/tools/" the `AccumulationRegions` class applies some postprocessing to the transport matrix simulation results to determine the boundary of the simulated accumulation regions.

Under "tm/plotting/" there are functions to: create an animation of the transport matrix simulation results, plot timeseries of the amount of simulated tracer contained in each of the accumulation regions, and to plot a map showing the simulated tracer concentrations at a requested time.

Some examples of how to use these tools are shown in the "main.m" script.

## Drifter tools
Under the "drifters/" folder, there are several tools available to analyse and plot the drifter data directly.
