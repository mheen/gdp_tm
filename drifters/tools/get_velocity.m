function [u_mean,v_mean,vel_mean,eke,grid] = get_velocity(gdp_data,drogued_status,grid_resolution)
% Gets ensemble mean drifter velocities and associated eddy kinetic energy.
%
% Input arguments:
% - gdp_data: data structure created by read_gdp_data
% - drogued_status [string]: "all", "drogued" or "undrogued"
% - grid_resolution [degrees]: resolution to construct a grid with
%
% Output:
% - u_mean: ensemble mean drifter velocity in the u-direction
% - v_mean: ensemble mean drifter velocity in the v-direction
% - vel_mean: ensemble mean absolute drifter velocity
% - eke: mean eddy kinetic energy
%
grid = GlobalGrid(grid_resolution);
[u,v,vel] = combine_velocities_in_each_grid_cell(gdp_data,drogued_status,grid);
[u_mean_1d,v_mean_1d,vel_mean_1d] = get_mean_velocities(u,v,vel,grid);
eke_1d = get_eke(u,v,u_mean_1d,v_mean_1d,grid);
u_mean = vec2mat(u_mean_1d,grid.lon_size);
v_mean = vec2mat(v_mean_1d,grid.lon_size);
vel_mean = vec2mat(vel_mean_1d,grid.lon_size);
eke = vec2mat(eke_1d,grid.lon_size);
end

function eke_mean = get_eke(u_1d,v_1d,u_mean,v_mean,grid)
eke_mean = nan(grid.lon_size*grid.lat_size,1);
for i = 1:(grid.lon_size*grid.lat_size)    
    ke = 0.5*((u_1d{i}-u_mean(i)).^2+(v_1d{i}-v_mean(i)).^2);
    if ~isempty(ke)
        eke_mean(i) = nansum(ke)/length(ke);
    end
end
end

function [u_mean,v_mean,vel_mean] = get_mean_velocities(u_1d,v_1d,vel_1d,grid)
u_mean = nan(grid.lon_size*grid.lat_size,1);
v_mean = nan(grid.lon_size*grid.lat_size,1);
vel_mean = nan(grid.lon_size*grid.lat_size,1);
for i = 1:(grid.lon_size*grid.lat_size)
    u_mean(i) = nanmean(u_1d{i});
    v_mean(i) = nanmean(v_1d{i});
    vel_mean(i) = nanmean(vel_1d{i});
end
end

function [u_1d,v_1d,vel_1d] = combine_velocities_in_each_grid_cell(gdp_data,drogued_status,grid)
% initialise velocity cells
u_1d = cell(grid.lon_size*grid.lat_size,1);
v_1d = cell(grid.lon_size*grid.lat_size,1);
vel_1d = cell(grid.lon_size*grid.lat_size,1);
% get velocities in each grid cell
drifter_ids = fieldnames(gdp_data);
for i = 1:length(drifter_ids)
    if isfield(gdp_data.(drifter_ids{i}),drogued_status)
        lon = gdp_data.(drifter_ids{i}).(drogued_status)(:,3);
        lat = gdp_data.(drifter_ids{i}).(drogued_status)(:,2);
        [lon_index,lat_index] = grid.get_index(lon,lat);
        index_1d = sub2ind([grid.lon_size,grid.lat_size],lon_index,lat_index);
        u = gdp_data.(drifter_ids{i}).(drogued_status)(:,4);
        v = gdp_data.(drifter_ids{i}).(drogued_status)(:,5);
        vel = gdp_data.(drifter_ids{i}).(drogued_status)(:,6);        
        for j = 1:length(lon)
            if lon_index(j) <= grid.lon_size && lat_index(j) <= grid.lat_size                
                u_1d{index_1d(j)} = [u_1d{index_1d(j)} u(j)];
                v_1d{index_1d(j)} = [v_1d{index_1d(j)} v(j)];
                vel_1d{index_1d(j)} = [vel_1d{index_1d(j)} vel(j)];
            end
        end
    end
end
end

