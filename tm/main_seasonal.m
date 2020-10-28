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
drogued_status = 'drogued';

config1 = ModelConfiguration('drogued_status',drogued_status,...
    'specify_time.type','month','specify_time.times','[1,2]');
config2 = ModelConfiguration('drogued_status',drogued_status,...
    'specify_time.type','month','specify_time.times','[3,4]');
config3 = ModelConfiguration('drogued_status',drogued_status,...
    'specify_time.type','month','specify_time.times','[5,6]');
config4 = ModelConfiguration('drogued_status',drogued_status,...
    'specify_time.type','month','specify_time.times','[7,8]');
config5 = ModelConfiguration('drogued_status',drogued_status,...
    'specify_time.type','month','specify_time.times','[9,10]');
config6 = ModelConfiguration('drogued_status',drogued_status,...
    'specify_time.type','month','specify_time.times','[11,12]');

tm1 = TransportMatrix(gdp_data,config1);
tm2 = TransportMatrix(gdp_data,config2);
tm3 = TransportMatrix(gdp_data,config3);
tm4 = TransportMatrix(gdp_data,config4);
tm5 = TransportMatrix(gdp_data,config5);
tm6 = TransportMatrix(gdp_data,config6);
tm = {tm1,tm2,tm3,tm4,tm5,tm6};

simulation = TmModel(config1,tm,'n_iterations_per_tm',1);

%% -------------------------------------------------------------------------
% Processing and plotting
% -------------------------------------------------------------------------
tracer_threshold = 2; % minimum tracer concentration in accumulation region
area_threshold = 50; % minimum number of grid cells in accumulation region
accumulation = AccumulationRegions(simulation,config1,tracer_threshold,area_threshold);

plot_year = 10;
plot_tm_simulation(config1,plot_year,1);

figure; plot_temporal_accumulations(config1,tracer_threshold,area_threshold);
