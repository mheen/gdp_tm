function animate_tm_simulation(config,requested_years,lon_range,lat_range,output_path)
% Creates an .avi animation of simulated tracer concentrations for a
% range of requested simulation years, for a specific longitude and
% latitude range.
%
% Input arguments:
% - config: ModelConfiguration object
% - requested_years: array of simulation years for which simulation tracer
%   concentration is plotted. units are years after simulation started. If
%   the requested year is not available in the simulation, the simulation
%   time nearest to it is plotted.
% - lon_range: array containing minimum and maximum longitude values
% - lat_range: array containing minimum and maximum latitude values
% - output_path: location where video is saved
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

% initialise movie
mov = VideoWriter(output_path);
mov.FrameRate = 2;
open(mov);
% open figure
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1,1,1])
colormaps;
c_ticks = [0:1:15];
for t = 1:length(requested_years)
    % get tracer at time t
    [time_index,passed_time] = TmModel.get_time_index(simulation,requested_years(t));
    tracer = simulation.tracer(:,:,time_index);
    tracer(tracer<1) = NaN;
    % create plot
    clf; hold on;
    if exist('m_map')
        m_proj('mercator','lat',lat_range,'lon',lon_range);
        h = m_pcolor(lon,lat,tracer);
        m_coast('patch',[0,0,0]);
        m_grid('xtick',[],'ytick',[]);
    else
        pcolor(lon,lat,tracer); shading flat;
        xlim(lon_range)
        ylim(lat_range)
    end
    if strcmpi(config.drogued_status,'drogued')
        title([num2str(passed_time),' years simulation: 15 m depth'],'fontsize',20,'fontname','times')
    elseif strcmpi(config.drogued_status,'undrogued')
        title([num2str(passed_time),' years simulation: ocean surface'],'fontsize',20,'fontname','times')
    end        
    % colorbar
    c = colorbar;
    caxis([c_ticks(1),c_ticks(end)]);
    c.Ticks = [0,1,2,5,10,15];
    colormap(c_map_tracer);
    c.Label.String = ['Tracer concentration [amount per ',num2str(config.dx),char(176),' cells]'];
    c.Label.FontSize = 18;
    c.Label.FontName = 'Times';
    c.FontSize = 18;
    c.FontName = 'Times';
    c.Location = 'SouthOutside';
    % logos: temp
    uwa = imread('E:\documents\04_writing_presentations\figures\logos\uwa.png');
    axes('pos',[0.78,0.40,0.07,0.07])
    l1 = imshow(uwa);
    uu = imread('E:\documents\04_writing_presentations\figures\logos\uu.png');
    axes('pos',[0.78,0.33,0.07,0.07])
    l2 = imshow(uu);
    % save frame to movie
    ff = getframe(gcf);
    writeVideo(mov,ff);
end
close(mov);    
end
