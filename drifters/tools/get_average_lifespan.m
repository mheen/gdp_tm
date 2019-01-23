function av_lifespan = get_average_lifespan(gdp_data,drogued_status)
% Calculates the average lifespan of all (dead) drifters in gdp_data with a
% specified drogued_status. Drifters that are still transmitting are not
% included in the average lifespan.
%
% Input arguments:
% - gdp_data: data structure as constructed by "read_gdp_data"
% - drogued_status [string]: string that can be: "all", "drogued", "undrogued"
%
% Output:
% - av_lifespan: average drifter lifespan [years]
drifter_ids = fieldnames(gdp_data);
death_codes = get_death_code(drifter_ids);
n = 1;
for i = 1:length(drifter_ids)
    if isfield(gdp_data.(drifter_ids{i}),drogued_status)
        if ~isempty(gdp_data.(drifter_ids{i}).(drogued_status))            
            if death_codes(i) ~= 0 % drifter is still transmitting
                start_time = gdp_data.(drifter_ids{i}).(drogued_status)(1,1);
                end_time = gdp_data.(drifter_ids{i}).(drogued_status)(end,1);
                lifespan(n) = end_time-start_time;
                n = n+1;
            end
        end
    end
end
av_lifespan = nanmean(lifespan)/365;
end
