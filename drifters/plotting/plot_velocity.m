function plot_velocity(u,v,vel,global_grid)
% Plots a filled contour map of the mean drifter velocities with the
% streamlines of drifter velocities.
%
% Input arguments:
% - u: mean u component of drifter velocity calculated by get_velocity
% - v: mean v component of drifter velocity calculated by get_velocity
% - vel: mean absolute drifter velocity calculated by get_velocity
% - grid: GlobalGrid object used to calculate mean drifter velocities
%
% Requirements:
% This function uses the m_map package by Rich Pawlowicz, which can be
% downloaded from:
% https://www.eoas.ubc.ca/~rich/map.html
% If you do not have this installed, a simple plot without coastlines is
% made.
%
colormaps;
streamline_density = 20;

% use only latitude values [-75,75]
% this needs to be done in advance to avoid NaNs in the streamline plot
l_lat = global_grid.lat >= -75 & global_grid.lat <= 75;
lat_size = sum(l_lat);
lon = repmat(global_grid.lon,lat_size,1);
lat = repmat(global_grid.lat(l_lat),global_grid.lon_size,1)';
u = u(l_lat,:);
v = v(l_lat,:);
vel = vel(l_lat,:);


% plot
hold on
if exist('m_map')
    m_proj('mercator','lat',[-75,75],'lon',[-180,180]);
    [m_x,m_y] = m_ll2xy(lon,lat); % convert lon and lat to map projection (for streamline plot)
    h = m_pcolor(lon,lat,vel);
    h2 = streamslice(m_x,m_y,u,v,streamline_density);
    set(h2,'color','k')
    m_coast('patch',[0,0,0]);
    m_grid('box','fancy','tickdir','in','col',[0,0,0],'fontname','times','fontsize',12);
else
    pcolor(lon,lat,vel);
    shading flat;
    h2 = streamslice(lon,lat,u,v,streamline_density);
    set(h2,'color','k')
    xlim([-180,180])
    ylim([-75,75])
    xlabel(['Longitude [',char(176),']'])
    ylabel(['Latitude [',char(176),']'])
    grid on
    box on
end
c = colorbar;
colormap(c_map_vel);
caxis([0,100]);
c.Ticks = [0,5,10,20,30,40,50,75,100];
c.Label.String = 'Absolute mean velocity [cm/s]';
c.Label.FontSize = 14;
c.Label.FontName = 'times';
set(gca,'fontname','times','fontsize',12);
set(gcf,'color',[1,1,1])
