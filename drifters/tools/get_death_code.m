function death_code = get_death_code(drifter_id)
% Returns the death code of a drifter.
%
% Input argument:
% - drifter_id: ID of the drifter (first structure of gdp_data)
%
% Output:
% - death_code: integer with death code, with the following meaning:
%   0: still alive
%   1: ran aground (beached)
%   2: picked up by vessel
%   3: stopped transmitting
%   4: sporadic transmissions
%   5: bad batteries
%   6: inactive status
death_codes = read_death_codes();
if ~iscell(drifter_id)
    drifter_id = {drifter_id};
end
for i = 1:length(drifter_id)
    id = str2num(drifter_id{i}(3:end));
    l_drifter = death_codes(:,1)==id;
    death_code(i) = death_codes(l_drifter,2);
end
end

function death_codes = read_death_codes()
json_file = fileread('dirs.json');
gdp_dirs = parse_json(json_file);
input_dir = gdp_dirs.dirs.input_dir;
meta_data_fields = gdp_dirs.dirs.meta_data_filenames;
for i = 1:length(meta_data_fields)
    meta_data_input_paths{i} = [input_dir,meta_data_fields{i}];
end
death_codes = [];
for i = 1:length(meta_data_input_paths)
    death_codes_single = read_meta(meta_data_input_paths{i});
    death_codes = [death_codes;death_codes_single];
end
end

function death_codes = read_meta(filename)
fid = fopen(filename,'r');
format_spec = '%d %*d %*d %*s %*s %*s %*f %*f %*s %*s %*f %*f %*s %*s %d';
raw = fscanf(fid,format_spec);
ids = raw([1:2:length(raw)-1]);
codes = raw([2:2:length(raw)]);
death_codes = [ids,codes];
fclose(fid);
end