function plot_tm_simulation(config,requested_year,l_accumulation)
% Plots simulated tracer concentration for a requested simulation year.
%
% Input arguments:
% - simulation: data structure containing simulation result, created by
%   running a transport matrix simulation using "TmModel"
% - requested_year: simulation year for which simulation tracer
%   concentration is plotted. units are years after simulation started. if
%   the requested year is not available in the simulation, the simulation
%   time nearest to it is plotted.
% - l_accumulation [logical]: flag that indicates if outlines of
%   accumulation regions should be plotted as well
%
% Requirements:
% This function uses the m_map package by Rich Pawlowicz, which can be
% downloaded from:
% https://www.eoas.ubc.ca/~rich/map.html
% If you do not have this installed, a simple plot without coastlines is
% made.
%

simulation = TmModel.load(config);

lon = repmat(simulation.grid.lon,length(simulation.grid.lat),1);
lat = repmat(simulation.grid.lat,length(simulation.grid.lon),1)';

[time_index,passed_time] = TmModel.get_time_index(simulation,requested_year);
tracer = simulation.tracer(:,:,time_index);
tracer(tracer<1) = NaN;

% reorganise so that Indian Ocean is plotted in middle of figure
center_lon = 98;
i_start = find(simulation.grid.lon==-180+center_lon)+1;
i_middle = length(simulation.grid.lon);
i_end = i_start-1;
ind=[i_start:i_middle,1:i_end]; % move right side to left
lon = lon(:,ind);
lon(lon<=-180+center_lon) = lon(lon<=-180+center_lon)+360;
lat = lat(:,ind);
tracer = tracer(:,ind);

colormaps;
c_ticks = [0:1:15];
hold on
if exist('m_map')
    m_proj('mercator','lat',[-75,75],'lon',center_lon);
    h = m_pcolor(lon,lat,tracer);
    m_coast('patch',[0,0,0]);
    m_grid('box','fancy','tickdir','in','col',[0,0,0]);
else
    pcolor(lon,lat,tracer); shading flat;
    xlim([-180,180])
    ylim([-75,75])
    xlabel(['Longitude [',char(176),']']);
    ylabel(['Latitude [',char(176),']'])
end
c = colorbar;
caxis([c_ticks(1),c_ticks(end)]);
c.Ticks = [0,1,2,5,10,15];
colormap(c_map_tracer);
c.Label.String = ['Tracer concentration [amount per ',num2str(config.dx),char(176),' cells]'];
c.Label.FontSize = 18;
c.Label.FontName = 'Times';
title([num2str(passed_time),' year simulation, ',config.drogued_status],...
    'fontsize',20,'fontname','times','fontweight','normal')
set(gcf,'color',[1,1,1])

if l_accumulation
    accumulation = AccumulationRegions.load(config,2,50);
    for i = 1:length(accumulation.names)
        if isfield(accumulation.outlines.(accumulation.names{i}){time_index},'lon')
            for j = 1:length(accumulation.outlines.(accumulation.names{i}){time_index}.lon)
                lon = accumulation.outlines.(accumulation.names{i}){time_index}.lon{j};
                lat = accumulation.outlines.(accumulation.names{i}){time_index}.lat{j};
                lon(lon<=-180+center_lon) = lon(lon<=-180+center_lon)+360;
                p = m_plot(lon,lat,'-k','linewidth',2);
            end
        end
    end
    l = legend(p,'Accumulation regions','location','northeast');
    l.FontSize = 12;
    l.FontName = 'times';
    
end
