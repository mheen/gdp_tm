function gdp_data_specific = get_gdp_data_specific_time(gdp_data,time_type,times)
% Returns GDP drifter data for specified month(s) or year(s) only.
%
% Input arguments:
% - gdp_data: data structure created by "read_gdp_data"
% - time_type [string]: "month" or "year"
% - times [array]: array of integers with months or years for which
%   GDP drifter data will be extracted
%
% Output:
% - gdp_data_month: data structure following the same format as gdp_data
%   (created by "read_gdp_data")
%
drifter_ids = fieldnames(gdp_data);
drogued_status = {'all','drogued','undrogued'};

for s = 1:length(drogued_status)
    for d = 1:length(drifter_ids)
        if ~isempty(gdp_data.(drifter_ids{d}).(drogued_status{s}))
            drifter_time = gdp_data.(drifter_ids{d}).(drogued_status{s})(:,1);
            if strcmpi(time_type,'year')
                l_times = ismember(year(drifter_time),times);
            elseif strcmpi(time_type,'month')
                l_times = ismember(month(drifter_time),times);
            else
                error('Unknown time_type requested. Valid options are "month" and "year".')
            end
            if any(l_times)
                gdp_data_specific.(drifter_ids{d}).(drogued_status{s}) = gdp_data.(drifter_ids{d}).(drogued_status{s})(l_times,:);
            end
        end
    end
end
end
