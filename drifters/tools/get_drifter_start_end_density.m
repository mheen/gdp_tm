function [drifter_start,drifter_end,grid] = get_drifter_start_end_density(gdp_data,drogued_status,grid_resolution)
% Gets density of start and end locations of drifters in square grid cells.
%
% Input arguments:
% - gdp_data: data structure as constructed by "read_gdp_data"
% - drogued_status [string]: string that can be: "all", "drogued", "undrogued"
% - grid_resolution [degrees]: value that specifies the resolution of the
%   square grid cells in which drifter start and end locations are noted
%
% Output:
% - drifter_start [lon,lat]: matrix with drifter start location densities
% - drifter_end [lon,lat]: matrix with drifter end location densities
% - grid: GlobalGrid object
%
grid = GlobalGrid(grid_resolution);
drifter_start = zeros(grid.lon_size,grid.lat_size);
drifter_end = zeros(grid.lon_size,grid.lat_size);
% get drifter start and end densities
drifter_ids = fieldnames(gdp_data);
for i = 1:length(drifter_ids)
    if isfield(gdp_data.(drifter_ids{i}),drogued_status)
        if ~isempty(gdp_data.(drifter_ids{i}).(drogued_status))
            lon_start = gdp_data.(drifter_ids{i}).(drogued_status)(1,3);
            lon_end = gdp_data.(drifter_ids{i}).(drogued_status)(end,3);
            lat_start = gdp_data.(drifter_ids{i}).(drogued_status)(1,2);            
            lat_end = gdp_data.(drifter_ids{i}).(drogued_status)(end,2);
            [lon_index_start,lat_index_start] = grid.get_index(lon_start,lat_start);
            [lon_index_end,lat_index_end] = grid.get_index(lon_end,lat_end);
            if ~isnan(lon_index_start) && ~isnan(lat_index_start)
                drifter_start(lon_index_start,lat_index_start) = drifter_start(lon_index_start,lat_index_start)+1;
            end
            if ~isnan(lon_index_end) && ~isnan(lat_index_end)
                drifter_end(lon_index_end,lat_index_end) = drifter_end(lon_index_end,lat_index_end)+1;
            end
        end
    end
end
end
