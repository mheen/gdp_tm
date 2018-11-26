function plot_drifter_density(drifter_density,grid,title_text)
% Plots drifter density on a log-scale.
%
% Input arguments:
% - drifter_density: matrix with drifter densities, created by the
%   "get_drifter_density" function
% - grid: GlobalGrid object, returned as output by the
%   "get_drifter_density" function
% - title_text: string, text that is displayed in the figure title
%
% Requirements:
% This function uses the m_map package by Rich Pawlowicz, which can be
% downloaded from:
% https://www.eoas.ubc.ca/~rich/map.html
% If you do not have this installed, a simple plot without coastlines is
% made.
%
colormaps;
lon = repmat(grid.lon,grid.lat_size,1);
lat = repmat(grid.lat,grid.lon_size,1)';

l_filter = drifter_density==0;
drifter_density(l_filter) = NaN;

% plot
hold on

if exist('m_map')
    m_proj('mercator','lat',[-75,75],'lon',0);
    h = m_pcolor(lon,lat,log10(drifter_density)');
    m_coast('patch',[0,0,0]);
    m_grid('box','fancy','tickdir','in','col',[0,0,0]);
else
    pcolor(lon,lat,log10(drifter_density)');
    shading flat;
    xlim([-180,180])
    ylim([-75,75])
    xlabel(['Longitude [',char(176),']'])
    ylabel(['Latitude [',char(176),']'])
end
c = colorbar;
caxis([log10(1),log10(10^4)]);
c_ticks = [1,10,10^2,10^3,10^4];
c.Ticks = log10(c_ticks);
c.TickLabels = {'1','10','10^2','10^3','10^4'};
c.FontSize = 10;
c.Label.String = ['Drifter location fixes [# per ',num2str(grid.resolution),char(176),' cells]'];
c.Label.FontSize = 12;
colormap(c_map_density);
set(gcf,'color',[1 1 1])
title(title_text,'fontsize',14);
