% -------------------------------------------------------------------------
% Transport matrix simulation based on GDP drifter data
% -------------------------------------------------------------------------
% This script runs tracer simulations using a transport matrix based on
% observed drifter locations from NOAA's Global Drifter Program (GDP).
%
% To run this script:
% 1. Load "gdp_data.mat". If you do not have this file yet, create it using
%     the function "../drifters/read_gdp_data". (This takes a while,
%     so you preferably only want to do this once.)
% 2. Make sure there is a "config.ini" file present in this folder. Adjust
%    the transport matrix and simulation settings in the .ini file as
%    desired.
% You are now ready to run this script!
%

% -------------------------------------------------------------------------
% Simulation
% -------------------------------------------------------------------------
config = ModelConfiguration('config.ini');

tm = TransportMatrix(gdp_data,config);

simulation = TmModel(config,tm);

% -------------------------------------------------------------------------
% Processing and plotting
% -------------------------------------------------------------------------
tracer_threshold = 2; % minimum tracer concentration in accumulation region
area_threshold = 50; % minimum number of grid cells in accumulation region
accumulation = AccumulationRegions(simulation,config,tracer_threshold,area_threshold);

plot_year = 10;
plot_tm_simulation(config,plot_year,1);

plot_temporal_accumulations(config,tracer_threshold,area_threshold);