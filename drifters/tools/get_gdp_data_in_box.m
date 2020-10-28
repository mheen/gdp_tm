function gdp_data_box = get_gdp_data_in_box(gdp_data,box)
% Finds drifters that have been within a specified boxed area
%
% Input arguments:
%  - gdp_data: data structure as constructed by "read_gdp_data"
%  - box: structure containing a field "lon" and a field "lat",
%         which are vectors containing the longitude and latitude
%         ranges of the boxed area
%
% Output:
%  - gdp_data_box: data structure containing data of drifters that have
%                  been in the box. structure is similar to the gdp_data
%                  structure, with data further separated as:
%                  "all": all drifter data in the box
%                  "before": drifter data before entering the box
%                  "during": drifter data while inside and re-entering box
%                  "after": drifter data after leaving the box
%
drogued_status = {'all','drogued','undrogued'};
drifter_ids = fieldnames(gdp_data);
for i = 1:length(drogued_status)
    for j = 1:length(drifter_ids)
        if isfield(gdp_data.(drifter_ids{j}),drogued_status{i})            
            data = gdp_data.(drifter_ids{j}).(drogued_status{i});
            drifter_lon = data(:,3);
            drifter_lat = data(:,2);
            if box.lon(2) < box.lon(1)
                in_lon_box_left = box.lon(2) < drifter_lon & drifter_lon < 180;
                in_lon_box_right = -180 < drifter_lon & drifter_lon < box.lon(1);
                in_lon_box = or(in_lon_box_left,in_lon_box_right);
            else
                in_lon_box = box.lon(1) < drifter_lon & drifter_lon < box.lon(2);
            end
            in_lat_box = box.lat(1) < drifter_lat & drifter_lat < box.lat(2);
            in_box = in_lon_box & in_lat_box;
            if any(in_box)
                gdp_data_box.(drifter_ids{j}).(drogued_status{i}).all = data;
                index_box = find(in_box);
                index_entry = index_box(1);
                index_exit = index_box(end);
                gdp_data_box.(drifter_ids{j}).(drogued_status{i}).before = data(1:index_entry-1,:);
                gdp_data_box.(drifter_ids{j}).(drogued_status{i}).during = data(index_entry:index_exit,:);
                gdp_data_box.(drifter_ids{j}).(drogued_status{i}).after = data(index_exit+1:end,:);
            end
        end
    end
end

if ~exist('gdp_data_box')
    gdp_data_box = NaN;
end
