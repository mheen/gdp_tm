classdef TransportMatrix < handle
    properties
        config % ModelConfiguration object
        grid % GlobalGrid object containing information on transport matrix grid
        tm % transport matrix
        tmn % normalised transport matrix (should be used for simulations)
        N % sum of each column in tm (drifters leaving a location)
        M % sum of each row in tm (drifters entering a location)        
        info % additional information regarding tm construction
    end
    
    methods
        function obj = TransportMatrix(gdp_data,config)
            % Get transport matrix from GDP drifter data. If a transport
            % matrix with the requested specifications (drogued status,
            % grid cell size dx and separation time dt) already exists in a
            % .mat file, this transport matrix is loaded. Otherwise, a new
            % transport matrix is constructed and saved to a .mat file.
            %
            % A non-normalised transport matrix is constructed first. Next,
            % sink locations are iteratively removed from the transport
            % matrix ("tm"). Finally, the transport matrix is normalised
            % ("tmn"). This normalised transport matrix should be used for
            % simulations.
            %
            % The transport matrix T is an n x n matrix, where n equals the
            % total number of grid cells. The columns j of T represent the
            % grid cell a drifter is in at time t and the rows i of T
            % represent the grid cell a drifter is in a time dt later.
            %
            % Additional data regarding the construction of the transport
            % matrix is stored in "info".
            % Information in "info" includes:
            % - sink locations
            % - locations of undefined grid cells (no drifter data)
            % - number of drifters used for construction
            % - number of drifter locations used for construction
            % - number of drifter locations in each month
            % - number of drifter locations in each year
            % 
            % Input arguments:
            % - gdp_data: data structure as constructed by "read_gdp_data"
            % - config: ModelConfiguration object
            %            
            try
                tm = obj.load(config);
                obj.config = tm.config;
                obj.grid = tm.grid;
                obj.tm = tm.tm;
                obj.tmn = tm.tmn;
                obj.N = tm.N;
                obj.M = tm.M;
                obj.info = tm.info;
            catch ME
                if(strcmp(ME.identifier,'TransportMatrix:fileNotFound'))
                    obj.config = config;
                    obj.grid = GlobalGrid(obj.config.dx);
                    % construct tm
                    if strcmpi(obj.config.specify_time.type,'all')
                        obj.get_tm(gdp_data)
                    else
                        obj.get_tm_specific_times(gdp_data)
                    end
                    obj.handle_sinks()
                    obj.get_normalised_tm()
                    obj.save()
                end
            end
        end
        
        function get_normalised_tm(obj)
            % Normalises the transport matrix by dividing the columns j of
            % the transport matrix T by the sum of each column N. This
            % ensures that the number of drifters is conserved after
            % multiplication with T.
            %
            obj.tmn = obj.tm;            
            obj.N = sum(obj.tmn,1);
            obj.M = sum(obj.tmn,2)';
            
            % normalise columns so that number of drifters is conserved
            index_normalise = find(obj.N~=0);
            for i = 1:length(index_normalise)
                j = index_normalise(i);
                obj.tmn(:,j) = obj.tmn(:,j)/obj.N(j);
            end
            
            % add undefined locations (no drifter data) to info
            obj.info.undefined.index = find(obj.N==0 & obj.M==0);
            [obj.info.undefined.j,obj.info.undefined.i] = ind2sub([obj.grid.lon_size,obj.grid.lat_size],obj.info.undefined.index);
            obj.info.undefined.lon = obj.grid.lon(obj.info.undefined.j);
            obj.info.undefined.lat = obj.grid.lat(obj.info.undefined.i);
        end
        
        function handle_sinks(obj)
            % Handles sink locations in the transport matrix T.
            % Sink locations are found according as explained below.
            %
            % The sum of the i-th row of T is M(i). M(i) is the total
            % number of drifter locations that enter grid cell i. The sum of
            % the j-th column of T is N(j). N(j) is the total number of
            % drifter locations that leave grid cell j.
            % In the simplest case, a sink is a grid cell k where M(k)~=0
            % and N(k)==0, indicating that drifters entered the cell
            % (M(k)~=0) but never left (N(k)==0).
            % Drifters can also have several observed locations in a sink
            % grid cell k. In that case, N(k)~=0. But in this case, the only
            % drifter locations that 'leave' grid cell k, enter back into it.
            % Then for these locations N(k)==T(k,k). This condition is also
            % satisfied when N(k)=0.
            % So a grid cell k is a sink location if it satisfies:
            % N(k)==T(k,k) & M(k)~=0,
            % the condition M(k)~=0 is added to exclude locations that
            % contain no drifter locations at all.
            %
            % Sinks in T are handled by one of two options (specified by
            % either 1 or 0 in the "sinks" parameter).
            %
            % 0: Sink locations are removed from T iteratively.
            % Sink locations k are removed by setting the entire row in the
            % transport matrix to zero: T(k,:)=0. This prevents drifters
            % from entering into sink locations at all.
            % Because removing a sink location sometimes only shifts the
            % sink, sinks are removed iteratively.
            %
            % 1: 'Fixing' sinks so that tracer remains in the sink
            % locations forever.
            % This is done by setting T(k,k) = 1, where k is a sink grid
            % cell.
            %
            obj.N = sum(obj.tm,1);
            obj.M = sum(obj.tm,2)';
            l_sink = obj.N == diag(obj.tm)' & obj.M~=0;
            sink_index = find(l_sink);            
            if obj.config.sinks == 0
                % iteratively remove sinks
                obj.info.sinks.n = 0;
                obj.info.sinks.iterations = 0;
                obj.info.sinks.index = [];
                while ~isempty(sink_index)
                    % get info
                    obj.info.sinks.n = obj.info.sinks.n+length(sink_index);
                    obj.info.sinks.iterations = obj.info.sinks.iterations+1;
                    obj.info.sinks.index = [obj.info.sinks.index,sink_index];
                    % remove sinks from tm
                    obj.tm(:,sink_index) = 0;
                    obj.tm(sink_index,:) = 0;
                    obj.N = sum(obj.tm,1);
                    obj.M = sum(obj.tm,2)';
                    l_sink = obj.N == diag(obj.tm)' & obj.M~=0;
                    sink_index = find(l_sink);
                end
            elseif obj.config.sinks == 1
                % fix sinks
                obj.info.sinks.n = length(sink_index);
                obj.info.sinks.index = sink_index;
                obj.tm(sink_index,sink_index) = 1;
            end
            % get sink locations for info
            [sink_j,sink_i] = ind2sub([obj.grid.lon_size,obj.grid.lat_size],obj.info.sinks.index);
            obj.info.sinks.lon = obj.grid.lon(sink_j);
            obj.info.sinks.lat = obj.grid.lat(sink_i);
        end
        
        function get_tm_specific_times(obj,gdp_data)
            % Creates a transport matrix T from gdp_data for specific
            % months or years. For each drifter, observed locations at a
            % time t that falls in the specified months or years and at a
            % time t+dt are noted. The grid cell that a drifter is in at t
            % is added to the corresponding column of T, the grid cell a
            % drifter is in at a time t+dt is added to the corresponding
            % row in T.
            %
            obj.check_gdp_data(gdp_data);
            
            grid_index_t0 = [];
            grid_index_t0dt = [];
            time_locations = [];
            n_drifters = 0;
            
            drifter_ids = fieldnames(gdp_data);
            for i = 1:length(drifter_ids)
                if ~isempty(gdp_data.(drifter_ids{i}).(obj.config.drogued_status))
                    drifter_time = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(:,1);
                    separation = obj.config.dt*4; % because drifter measurements are returned every 6 hours (4 times a day)
                    if strcmpi(obj.config.specify_time.type,'month')
                        times = str2num(obj.config.specify_time.times);
                        l_times = ismember(month(drifter_time),times);
                    elseif strcmpi(obj.config.specify_time.type,'year')
                        times = str2num(obj.config.specify_time.times);
                        l_times = ismember(year(drifter_time),times);
                    else
                        error('Unknown type of specify_time requested. Valid options are: "all", "month", "year".')
                    end
                    if any(l_times)
                        % get time indices of locations
                        time_index_t0 = find(l_times);
                        time_index_t0dt = time_index_t0+separation;
                        i_inrange = find(time_index_t0dt<=length(drifter_time)); % keep only times where t0+dt is in range of drifter measurement
                        if any(i_inrange)
                            n_drifters = n_drifters+1;
                            time_index_t0 = time_index_t0(i_inrange);
                            time_index_t0dt = time_index_t0dt(i_inrange);
                            % get lat and lon of locations
                            lon_t0 = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0,3);
                            lat_t0 = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0,2);
                            lon_t0dt = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0dt,3);
                            lat_t0dt = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0dt,2);
                            % get lat and lon indexes of locations
                            l_no_nan = ~isnan(lon_t0) & ~isnan(lat_t0) & ~isnan(lon_t0dt) & ~isnan(lat_t0dt);
                            [lon_index_t0,lat_index_t0] = obj.grid.get_index(lon_t0(l_no_nan),lat_t0(l_no_nan));
                            [lon_index_t0dt,lat_index_t0dt] = obj.grid.get_index(lon_t0dt(l_no_nan),lat_t0dt(l_no_nan));
                            % get grid index of locations
                            grid_index_t0_append = sub2ind([obj.grid.lon_size,obj.grid.lat_size],lon_index_t0,lat_index_t0);
                            grid_index_t0dt_append = sub2ind([obj.grid.lon_size,obj.grid.lat_size],lon_index_t0dt,lat_index_t0dt);
                            grid_index_t0 = [grid_index_t0;grid_index_t0_append];
                            grid_index_t0dt = [grid_index_t0dt;grid_index_t0dt_append];
                            % get time of drifter locations (for info)
                            time_locations_append = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0dt(l_no_nan),1);
                            time_locations = [time_locations;time_locations_append];
                        end
                    end
                end
            end
            locations = ones(length(grid_index_t0dt),1);
            obj.tm = sparse(grid_index_t0dt,grid_index_t0,locations,obj.grid.lon_size*obj.grid.lat_size,obj.grid.lon_size*obj.grid.lat_size);
            
            obj.get_info(n_drifters,time_locations);            
        end
        
        function get_tm(obj,gdp_data)
            % Creates a transport matrix T from gdp_data. For each drifter,
            % observed locations at a time t and a time t+dt are noted. The
            % grid cell that a drifter is in at t is added to the
            % corresponding column of T, the grid cell a drifter is in at a
            % time t+dt is added to the corresponding row in T.
            %
            % IMPORTANT: This function assumes that all drifter
            % measurements are returned every 6 hours. This condition is
            % checked before constructing the transport matrix.
            %
            % Input arguments:
            % - gdp_data: data structure as constructed by "read_gdp_data"
            %
            obj.check_gdp_data(gdp_data);
            
            grid_index_t0 = [];
            grid_index_t0dt = [];
            time_locations = [];
            n_drifters = 0;
            
            drifter_ids = fieldnames(gdp_data);
            for i = 1:length(drifter_ids)
                if ~isempty(gdp_data.(drifter_ids{i}).(obj.config.drogued_status))
                    n_drifters = n_drifters+1;
                    % get time indices of locations
                    separation = obj.config.dt*4; % because drifter measurements are returned every 6 hours (4 times a day)
                    time_index_t0 = [1:length(gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(:,1))];
                    time_index_t0dt = time_index_t0+separation;
                    i_cutoff = time_index_t0(end)-separation; % cut-off where i_t0dt will exceed available time;
                    time_index_t0 = time_index_t0(1:i_cutoff);
                    time_index_t0dt = time_index_t0dt(1:i_cutoff);
                    % get lat and lon of locations
                    lon_t0 = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0,3);
                    lat_t0 = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0,2);
                    lon_t0dt = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0dt,3);
                    lat_t0dt = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0dt,2);
                    % get lat and lon indexes of locations
                    l_no_nan = ~isnan(lon_t0) & ~isnan(lat_t0) & ~isnan(lon_t0dt) & ~isnan(lat_t0dt);                    
                    [lon_index_t0,lat_index_t0] = obj.grid.get_index(lon_t0(l_no_nan),lat_t0(l_no_nan));
                    [lon_index_t0dt,lat_index_t0dt] = obj.grid.get_index(lon_t0dt(l_no_nan),lat_t0dt(l_no_nan));
                    % get grid index of locations
                    grid_index_t0_append = sub2ind([obj.grid.lon_size,obj.grid.lat_size],lon_index_t0,lat_index_t0);
                    grid_index_t0dt_append = sub2ind([obj.grid.lon_size,obj.grid.lat_size],lon_index_t0dt,lat_index_t0dt);
                    grid_index_t0 = [grid_index_t0;grid_index_t0_append];
                    grid_index_t0dt = [grid_index_t0dt;grid_index_t0dt_append];
                    % get time of drifter locations (for info)
                    time_locations_append = gdp_data.(drifter_ids{i}).(obj.config.drogued_status)(time_index_t0dt(l_no_nan),1);
                    time_locations = [time_locations;time_locations_append];
                end
            end
            locations = ones(length(grid_index_t0dt),1);
            obj.tm = sparse(grid_index_t0dt,grid_index_t0,locations,obj.grid.lon_size*obj.grid.lat_size,obj.grid.lon_size*obj.grid.lat_size);
            
            obj.get_info(n_drifters,time_locations);
        end
        
        function get_info(obj,n_drifters,time_locations)
            % Gets some information on the construction of the transport
            % matrix T, to store in "info".
            %
            obj.info.n_drifters = n_drifters;
            obj.info.n_locations = length(time_locations);
            month_locations = month(time_locations);
            unique_months = unique(month_locations,'sorted');
            obj.info.months = unique_months;
            obj.info.n_months = histc(month_locations,unique_months);
            year_locations = year(time_locations);
            unique_years = unique(year_locations,'sorted');
            obj.info.years = unique_years;
            obj.info.n_years = histc(year_locations,unique_years);
        end
        
        function save(obj)
            % Writes TransportMatrix object to "tm" data structure and
            % saves in a .mat file.
            %
            [tm_dir,tm_file] = obj.get_path(obj.config);
            if ~exist(tm_dir,'dir')
                mkdir(tm_dir);
            end
            tm.config.drogued_status = obj.config.drogued_status;
            tm.config.dx = obj.config.dx;
            tm.config.dt = obj.config.dt;
            tm.config.sinks = obj.config.sinks;
            tm.grid = obj.grid;
            tm.tm = obj.tm;
            tm.tmn = obj.tmn;
            tm.N = obj.N;
            tm.M = obj.M;
            tm.info = obj.info;            
            save([tm_dir,tm_file],'tm','-v7.3');
        end
    end
    
    methods(Static)
        function [tm_dir,tm_file] = get_path(config)
            % Static function that returns the directory and file name
            % where a transport matrix is stored, based on specific
            % parameters of the transport matrix.
            %
            % Input arguments:
            % - config: ModelConfiguration object
            %
            % Output:
            % - tm_dir: directory where .mat file with transport matrix is
            %   stored
            % - tm_file: filename where .mat file with transport matrix is
            %   stored
            %
            dx_description = ['_dx',strrep(num2str(config.dx),'.','')];
            dt_description = ['_dt',num2str(config.dt)];
            sink_description = ['_s',num2str(config.sinks)];
            if strcmpi(config.specify_time.type,'all')
                time_description = '';
            elseif strcmpi(config.specify_time.type,'month') || strcmpi(config.specify_time.type,'year')
                times = str2num(config.specify_time.times);
                times_str = num2str(times(1));
                for t = 2:length(times)
                    times_str = [times_str,'-',num2str(times(t))];
                end
                time_description = ['_',config.specify_time.type(1),times_str];
            else
                error('Unknown specify_time type requested. Valid options are: "all", "month", "year".')
            end
            
            run_dir = fileparts(mfilename('fullpath'));
            tm_dir = [run_dir,'/output/',config.drogued_status,dx_description,dt_description,sink_description,time_description,'/'];
            tm_file = 'tm.mat';
        end
        
        function tm = load(config)
            [tm_dir,tm_file] = TransportMatrix.get_path(config);
            if exist([tm_dir,tm_file])
                load([tm_dir,tm_file]);
            else
                error('TransportMatrix:fileNotFound','Transport matrix does not exist yet. Please create it first.');
            end
        end
        
        function check_gdp_data(gdp_data)
            % Checks if GDP drifter data is uniformly spaced 6 hours apart.
            % Returns error message if this is not the case (essential in
            % constructing the transport matrix).
            %
            % Input arguments:
            % - gdp_data:  data structure as constructed by "read_gdp_data"
            %
            drifter_ids = fieldnames(gdp_data);
            for i = 1:length(drifter_ids)
                dt = gdp_data.(drifter_ids{i}).all(2:end,1)-gdp_data.(drifter_ids{i}).all(1:end-1,1);
                uniques_dt = length(unique(dt));
                if uniques_dt == 1
                    if dt ~= 0.25
                        error(['Fatal for transport matrix derivation: Data for drifter ',drifter_ids{i},' not spaced 6 hours apart.'])
                    end
                elseif uniques_dt > 1
                    error(['Fatal for transport matrix derivation: Non uniform spacing in time of data for drifter ',drifter_ids{i},'.'])
                end
            end
        end        
    end
end
