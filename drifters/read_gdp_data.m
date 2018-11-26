function gdp_data = read_gdp_data()
% Reads GDP drifter data from "buoydata_*.dat" and "dirfl_*.dat" files and
% saves in a "*.mat" file.
% These files are available through NOAA's GDP program. For information on
% how to download this data see:
% http://www.aoml.noaa.gov/phod/gdp/interpolated/data/all.php
% or download directly from the FTP-server:
% ftp://ftp.aoml.noaa.gov/pub/phod/buoydata/
% The file locations are read from "gdp/drifters/dirs.json", which should
% be structured as outlined in the README file. The output file path is
% also retrieved from "dirs.json".
% 
% The final "gdp_data.mat" file contains a data structure which is
% organised as follows:
% Level 1. "gdp_data"
%       Level 2. drifter ID (in the form of "IDxxxxxxxx")
%             Level 3. "raw" (matrix with full "buoydata_*.dat" content)
%             Level 3. "all" (matrix with full drifter data)
%             Level 3. "drogue_release_time" (MATLAB datenum)
%             Level 3. "drogued" (matrix with drogued drifter data)
%             Level 3. "undrogued" (matrix with undrogued drifter data)
% The matrices in level 3 "all", "drogued" and "undrogued" consist of 6
% columns representing:
% Column 1. time (MATLAB datenum)
% Column 2. latitude
% Column 3. longitude
% Column 4. u-velocity component
% Column 5. v-velocity component
% Column 6. absolute velocity
%
% Input arguments:
%  - "dirs.json" in "gdp/drifters" directory set up as detailed in README
% Output arguments:
% - gdp_data: data structure as described above
% - .mat file containing gdp_data saved to output path in "dirs.json"
%
addpath('../')
% load file input and output paths from "dirs.json"
[data_input_paths,meta_input_paths,output_path] = read_filenames_from_json('dirs.json');

gdp_data = read(data_input_paths,meta_input_paths);
gdp_data = replace_missing_values_with_nan(gdp_data);
save(output_path,'gdp_data','-v7.3')

function [data_input_paths,meta_input_paths,output_path] = read_filenames_from_json(input_path)
json_file = fileread(input_path);
gdp_dirs = parse_json(json_file);
input_dir = gdp_dirs.dirs.input_dir;
output_path = gdp_dirs.dirs.output_path;
data_fields = gdp_dirs.dirs.data_filenames;
meta_data_fields = gdp_dirs.dirs.meta_data_filenames;
if length(data_fields) ~= length(meta_data_fields)
    error('Number of files containing GDP data does not match number of files containing meta data.');
else
    for i = 1:length(data_fields)
        data_input_paths{i} = [input_dir,data_fields{i}];
        meta_input_paths{i} = [input_dir,meta_data_fields{i}];
    end
end

function gdp_data = read(data_input_paths,meta_input_paths)
raw_data = [];
raw_meta = [];
for i = 1:length(data_input_paths)
    raw_data_single = read_data(data_input_paths{i});
    raw_meta_data_single = read_meta(meta_input_paths{i});
    raw_data = [raw_data;raw_data_single];
    raw_meta = [raw_meta;raw_meta_data_single];
end
[meta_drifter_id,drogue_release_time] = get_drogue_release_time_from_meta(raw_meta);
gdp_data = separate_data_into_all_drogued_undrogued(raw_data,meta_drifter_id,drogue_release_time);

function raw_data = read_data(filename)
fid = fopen(filename,'r');
format_spec = '%d %d %f %d %f %f %f %f %f %f %f %f %f';
data_size = [13 Inf];
raw_data = fscanf(fid,format_spec,data_size);
raw_data = raw_data';
fclose(fid);

function raw_meta = read_meta(filename)
fid = fopen(filename,'r');
format_spec = '%d %*d %*d %*s %*s %*s %*f %*f %*s %*s %*f %*f %s %s %*d';
size_meta = [16 Inf];
raw_meta = fscanf(fid,format_spec,size_meta);
raw_meta = raw_meta';
fclose(fid);

function [meta_drifter_id,drogue_release_time] = get_drogue_release_time_from_meta(raw_meta)
% drifter IDs
meta_drifter_id = raw_meta(:,1);
% drogue release time
str_drogue_release_time = char(raw_meta(:,2:end));
drogue_release_time = datenum(str_drogue_release_time,'yyyy/mm/ddHH:MM'); % 0000/00/0000:00 converts to -31
l_not_lost = drogue_release_time == -31; % find drifters that have not yet lost their drogue
drogue_release_time(l_not_lost) = datenum(date,'dd-mm-yyyy')+datenum(20,0,0); % replace for date 20 years in future

function gdp_data = separate_data_into_all_drogued_undrogued(raw_data,meta_drifter_id,drogue_release_time)
drifter_ids = unique(raw_data(:,1));
for i = 1:length(drifter_ids)
    fieldname = ['ID',num2str(drifter_ids(i))];
    % sort data per drifter
    gdp_data.(fieldname).raw = raw_data(raw_data(:,1)==drifter_ids(i),:);
    % get time
    year = gdp_data.(fieldname).raw(:,4);
    month = gdp_data.(fieldname).raw(:,2);
    day = floor(gdp_data.(fieldname).raw(:,3));
    hour = (gdp_data.(fieldname).raw(:,3)-day)*24;
    minute = zeros(size(gdp_data.(fieldname).raw,1),1);
    second = zeros(size(gdp_data.(fieldname).raw,1),1);
    datevector = [year month day hour minute second];
    gdp_data.(fieldname).all(:,1) = datenum(datevector);
    % get position
    lat = gdp_data.(fieldname).raw(:,5);
    lon = gdp_data.(fieldname).raw(:,6);
    % convert longitude [0,360] to [-180,180]
    if any(lon > 180)
        i360 = lon > 180;
        lon(i360) = lon(i360)-360;
    end
    gdp_data.(fieldname).all(:,2:3) = [lat lon];
    % get velocity
    u = gdp_data.(fieldname).raw(:,8);
    v = gdp_data.(fieldname).raw(:,9);
    vel = gdp_data.(fieldname).raw(:,10);
    gdp_data.(fieldname).all(:,4:6) = [u v vel];    
    % match drogue release time to drifter
    gdp_data.(fieldname).drogue_release_time = drogue_release_time(drifter_ids(i)== meta_drifter_id);
    % separate data with and without drogue
    i_drogue = gdp_data.(fieldname).all(:,1) < gdp_data.(fieldname).drogue_release_time;
    i_undrogue = gdp_data.(fieldname).all(:,1) >= gdp_data.(fieldname).drogue_release_time;
    gdp_data.(fieldname).drogued = gdp_data.(fieldname).all(i_drogue,:);
    gdp_data.(fieldname).undrogued = gdp_data.(fieldname).all(i_undrogue,:);    
end

function gdp_data = replace_missing_values_with_nan(gdp_data)
% Replaces NOAA's missing drifter values with NaN values.
% Missing values that NOAA uses are:
% latitude: 999.999
% longitude: 639.999
% velocity: 999.999
%
drogued_status = {'all','drogued','undrogued'};
for j = 1:length(drogued_status)    
    drifter_ids = fieldnames(gdp_data);    
    for i = 1:length(drifter_ids)
        if ~isempty(gdp_data.(drifter_ids{i}).(drogued_status{j}))
            % replace longitude data
            lon = gdp_data.(drifter_ids{i}).(drogued_status{j})(:,3);
            lon(lon==639.999) = NaN;
            gdp_data.(drifter_ids{i}).(drogued_status{j})(:,3) = lon;
            % replace latitude data
            lat = gdp_data.(drifter_ids{i}).(drogued_status{j})(:,2);
            lat(lat==999.999) = NaN;
            gdp_data.(drifter_ids{i}).(drogued_status{j})(:,2) = lat;
            % replace velocity data
            u = gdp_data.(drifter_ids{i}).(drogued_status{j})(:,4);
            v = gdp_data.(drifter_ids{i}).(drogued_status{j})(:,5);
            vel = gdp_data.(drifter_ids{i}).(drogued_status{j})(:,6);
            u(u==999.999) = NaN;
            v(v==999.999) = NaN;
            vel(vel==999.999) = NaN;
            gdp_data.(drifter_ids{i}).(drogued_status{j})(:,4) = u;
            gdp_data.(drifter_ids{i}).(drogued_status{j})(:,5) = v;
            gdp_data.(drifter_ids{i}).(drogued_status{j})(:,6) = vel;
        end
    end
end
