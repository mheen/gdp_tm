function [drifter_density,grid] = get_drifter_density(gdp_data,drogued_status,grid_resolution)
% Determines global drifter densities in square grid cells.
%
% Input arguments:
% - gdp_data: data structure as constructed by "read_gdp_data"
% - drogued_status: string that can be: "all", "drogued" or "undrogued"
% - grid_resolution: (degrees) value that specifies the resolution of the
%   square grid cells in which drifter density is calculated.
%
% Output:
% - drifter_density[longitude,latitude]: matrix with drifter density
% - grid: GlobalGrid object
%
addpath('../')
grid = GlobalGrid(grid_resolution);
drifter_density = zeros(grid.lon_size,grid.lat_size);

% calculate drifter density
drifter_ids = fieldnames(gdp_data);
for i = 1:length(drifter_ids)
    if ~isempty(gdp_data.(drifter_ids{i}).(drogued_status))
        lon = gdp_data.(drifter_ids{i}).(drogued_status)(:,3);
        lat = gdp_data.(drifter_ids{i}).(drogued_status)(:,2);
        [lon_index,lat_index] = grid.get_index(lon,lat);
        for j = 1:size(gdp_data.(drifter_ids{i}).(drogued_status),1)
            if ~isnan(lon_index(j)) && ~isnan(lat_index(j))
                drifter_density(lon_index(j),lat_index(j)) = drifter_density(lon_index(j),lat_index(j))+1;
            end
        end
    end
end