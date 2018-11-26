function plot_temporal_accumulations(config,tracer_threshold,area_threshold)
% Plots temporal evolution of percentage of tracer contained in
% accumulation regions. AccumulationRegions object is loaded or created
% based on the simulation and specified tracer and area threshold.
%
% Input arguments:
% - simulation: data structure containing simulation result, created by
%   running a transport matrix simulation using "TmModel"
% - tracer_threshold: tracer amount per grid cell that is set as a minimum
%   value to determine accumulation regions
% - area_threshold: minimum number of adjacent grid cells that make up an
%   accumulation region
% Loads:
% - accumulation: data structure or object containing information about
%   accumulation regions, created by AccumulationRegions
%
accumulation = AccumulationRegions.load(config,tracer_threshold,area_threshold);
time = accumulation.time/config.days_in_year;

hold on
max_percentage = [];
for i = 1:length(accumulation.names)
    tracer_percentage = cell2mat(accumulation.tracer.(accumulation.names{i}))/accumulation.total_tracer*100;
    max_percentage = [max_percentage,max(tracer_percentage)];
    plot(time,tracer_percentage,'linewidth',2);
end
l = legend(accumulation.names,'location','northeast');
l.FontSize = 12;
l.FontName = 'times';
xlabel('Time [years]','fontsize',18,'fontname','times');
ylabel('Tracer [%]','fontsize',18,'fontname','times');
xlim([0,max(time)]);
ylim([0,ceil(max(max_percentage))]);
grid on
title(['Temporal evolution of accumulation regions: ',config.drogued_status],'fontsize',20,'fontname','times');
set(gcf,'color',[1,1,1])
end
