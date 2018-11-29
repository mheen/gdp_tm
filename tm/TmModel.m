classdef TmModel < handle
    properties
        config % ModelConfiguration object
        initial % initial condition
        grid % GlobalGrid object containing information on transport matrix grid        
        tm % TransportMatrix object
        tracer % simulated tracer concentration
        time % [days] simulation time        
    end
    
    methods
        function obj = TmModel(config,tm)
            % Transport matrix model. Runs transport matrix simulation with
            % settings as specified in "config" and using transport matrix
            % contained in "tm".
            %
            % Input arguments:
            % - config: ModelConfiguration object, contains model and
            %   simulation settings read from a "config.ini" file
            % - tm: TransportMatrix object, transport matrix created from
            %   GDP drifter data, with settings contained in config
            %
            % Simulation output:
            % - tracer [lat,lon,time]: 3d matrix containing simulated
            %   tracer concentration
            % - time [days]
            % - grid: GlobalGrid object containing lon and lat
            %
            try
                simulation = obj.load(config);
                obj.config = simulation.config;                
                obj.grid = simulation.grid;
                obj.tracer = simulation.tracer;
                obj.time = simulation.time;                
            catch ME
                if(strcmp(ME.identifier,'TmModel:fileNotFound'))
                    obj.config = config;
                    obj.tm = tm;
                    obj.grid = GlobalGrid(obj.config.dx);            
                    % simulation
                    obj.get_initial_condition()
                    obj.run()                        
                    % save simulation and config.ini
                    obj.save()
                    [simulation_dir,~] = obj.get_path(config);
                    config.write_config(simulation_dir)
                end
            end
        end
        
        function get_initial_condition(obj)
            % Creates initial condition, depending on requested settings.
            %
            if strcmpi(obj.config.initial.type,'global_uniform')
                obj.get_global_uniform();
            elseif strcmpi(obj.config.initial.type,'uniform')
                obj.get_uniform_range();
            elseif strcmpi(obj.config.initial.type,'jambeck2015')
                error('Jambeck et al. (2015) initial condition type not yet available.')
            elseif strcmpi(obj.config.initial.type,'lebreton2017')
                error('Lebreton et al. (2017) initial condition type not yet available.')
            end
        end
        
        function get_global_uniform(obj)
            % Creates a global uniform initial condition, where tracer is
            % spread evenly over the global oceans.
            %
            obj.initial = obj.config.initial.tracer_per_grid_cell*ones(1,obj.grid.lon_size*obj.grid.lat_size);
            obj.initial(obj.tm.info.undefined.index) = 0;
        end
        
        function get_uniform_range(obj)
            % Creates a uniform initial condition, where tracer is spread
            % evenly between a specified longitude and latitude range.
            %
            lon1d = str2double(obj.config.initial.lon);
            lat1d = str2double(obj.config.initial.lat);
            lon2d = repmat(lon1d,length(lat1d),1);
            lat2d = repmat(lat1d,length(lon1d),1)';
            lon = lon2d(:);
            lat = lat2d(:);
                        
            [index_lon,index_lat] = obj.grid.get_index(lon,lat);
            index = sub2ind([obj.grid.lon_size,obj.grid.lat_size],index_lon,index_lat);
            
            obj.initial = zeros(1,obj.grid.lon_size*obj.grid.lat_size);
            obj.initial(index) = 1*obj.config.initial.tracer_per_grid_cell;
            obj.initial(obj.tm.info.undefined.index) = 0;
        end
                
        function run(obj)
            % Runs the simulation.
            %
            n_iterations = obj.get_iterations(obj.config.run_time);
            n_total_output = floor((n_iterations-1)/obj.config.output_interval)+1;
            % initialise matrices
            obj.time = nan(n_total_output,1);
            obj.time(1) = 0;
            tracer_2d = nan(n_total_output,size(obj.tm.tmn,2));
            tracer_2d(1,:) = obj.initial;
            old_state = obj.initial';
            n_output = 2;
            n_source_iterations = obj.add_sources_during_simulation();
            % run
            for i = 2:n_iterations
                if 2<i && i<=n_source_iterations
                    new_state = obj.tm.tmn*(old_state+obj.initial');
                else
                    new_state = obj.tm.tmn*old_state;
                end
                if mod(i-1,obj.config.output_interval) == 0
                    tracer_2d(n_output,:) = new_state;
                    obj.time(n_output) = (i-1)*obj.config.dt;
                    n_output = n_output+1;
                end
                old_state = new_state;
            end
            obj.tracer = obj.convert_2d_to_3d(tracer_2d,obj.grid);
            
            if obj.config.normalise
                obj.normalise_tracer()
            end
        end
        
        function normalise_tracer(obj)
            % Normalises tracer amount in a grid cell c(x,y) by multiplying
            % with the total number of grid cells N, divided by the total
            % amount of tracer in the simulation C,
            % TAF(x,y) = c(x,y)*N/C.
            % TAF is the tracer amplification factor as defined in van
            % Sebille et al. (2012): "If, for example, in some grid box
            % TAF(x,y)=20, then twenty times more tracer is found within
            % that grid box than if all tracer is uniformly distributed
            % over the global ocean."
            %
            N = obj.grid.lon_size*obj.grid.lat_size;
            for t = 1:size(obj.tracer,3)
                C = sum(sum(obj.tracer(:,:,t)));
                TAF = obj.tracer(:,:,t)*N/C;
                obj.tracer(:,:,t) = TAF;
            end
        end
        
        function n_source_iterations = add_sources_during_simulation(obj)
            % Checks if sources should be added during the simulation and
            % for how many iterations.
            %
            % Output:
            % - n_source_iterations [int]: number of model run iterations
            %   for which sources should be added, if no sources should be
            %   added n_source_iterations = 0 is returned.
            %
            if strcmpi(obj.config.initial.type,'jambeck2015') || strcmpi(obj.config.initial.type,'lebreton2017')
                if obj.config.initial.years_to_add_sources
                    n_source_iterations = obj.get_iterations(obj.config.initial.years_to_add_sources);
                else
                    n_source_iterations = 0;
                end
            else
                n_source_iterations = 0;
            end
        end
        
        function n_model_iterations = get_iterations(obj,run_time)
            % Calculates the number of required simulation iterations based
            % on the timestep dt (in days), the number of days in a year
            % (needs to be a multiple of the timestep dt) and the requested
            % run time (in years).
            %
            % Output:
            % - n_model_iterations [int]
            %
            run_time_days = run_time*obj.config.days_in_year;
            if mod(run_time_days,obj.config.dt) ~=0
                error(['Model run_time is not a multiple of the timestep dt.',...
                    'Check if days_in_year in config is a multiple of dt.'])
            else
                n_model_iterations = run_time_days/obj.config.dt+1;
            end
        end
                       
        function defined_sources = move_source_locations(obj,sources)
           % Checks if source locations are located at an undefined
           % location in the transport matrix. If they are, checks if there
           % is a defined location available within 1 grid cell. If so,
           % moves the source location to that defined grid cell.
           %
           defined_sources = sources;
           i_source = find(sources~=0);
           for i = 1:length(i_source)
               if ~isempty(find(i_source(i)==obj.undefined.index,1))
                   % k is the index of a grid cell that has a source, but is
                   % undefined in the transport matrix
                   k = i_source(i);
                   % indexes of direct neighbours of k
                   i_above = k-1;                   
                   i_below = k+1;                   
                   i_left = k+obj.grid.lat_size;                   
                   i_right = k-obj.grid.lat_size;
                   k_neighbours = [i_above,i_below,i_left,i_right];                   
                   % check if direct neighbours of k are defined or not
                   def_above = isempty(find(i_above==obj.undefined.index,1));
                   def_below = isempty(find(i_below==obj.undefined.index,1));
                   def_left = isempty(find(i_left==obj.undefined.index,1));
                   def_right = isempty(find(i_right==obj.undefined.index,1));
                   def_k_neighbours = [def_above,def_below,def_left,def_right];
                   % determine how to move k
                   if ~def_above && def_below && def_left && def_right
                       % assume coastline is above, move k below
                       defined_sources(k) = 0;
                       defined_sources(i_below) = sources(k);
                   elseif ~def_below && def_above && def_left && def_right
                       % assume coastline is below, move k above
                       defined_sources(k) = 0;
                       defined_sources(i_above) = sources(k);
                   elseif ~def_left && def_above && def_below && def_right
                       % assume coastline is left, move k right
                       defined_sources(k) = 0;
                       defined_sources(i_right) = sources(k);
                   elseif ~def_right && def_above && def_below && def_left
                       % assume coastline is right, move k left
                       defined_sources(k) = 0;
                       defined_sources(i_left) = sources(k);
                   elseif sum(def_k_neighbours)==1
                       % move k to only defined direct neighbour
                       def_nb = find(def_k_neighbours);
                       defined_sources(k) = 0;
                       defined_sources(k_neighbours(def_nb)) = sources(k);
                   elseif sum(def_k_neighbours)==2 || sum(def_k_neighbours)==4
                       % determine direct neighbours with max number of
                       % defined neighbours, if equal choose randomly
                       def_nbs = find(def_k_neighbours);
                       n_def_next_neighbours = zeros(1,length(def_nbs));
                       for j = 1:length(def_nbs)
                          m = k_neighbours(def_nbs(j)); % defined direct neighbour of k
                          next_neighbours = [m-1,m+1,...
                              m-obj.grid.lat_size,m-obj.grid.lat_size-1,m-obj.grid.lat_size+1,...
                              m+obj.grid.lat_size,m+obj.grid.lat_size-1,m+obj.grid.lat_size+1];
                          n_def = 0;
                          for l = 1:length(next_neighbours)
                              if isempty(find(next_neighbours(l)==obj.undefined.index,1))
                                  n_def = n_def+1;
                              end
                          end
                          n_def_next_neighbours(j) = n_def;
                       end
                       j_max = find(n_def_next_neighbours==max(n_def_next_neighbours));
                       if length(j_max)==1
                           k_new = k_neighbours(def_nbs(j_max));
                       else
                           j_random = j_max(randi(length(j_max)));
                           k_new = k_neighbours(def_nbs(j_random));
                       end
                       defined_sources(k) = 0;
                       defined_sources(k_new) = sources(k);
                   end                   
               end                   
           end
        end
        
        function save(obj)
            % Writes the simulation results to a data structure and saves
            % to a .mat file. The .mat file is saved to the directory that
            % also contains the relevant transport matrix.
            %
            [simulation_dir,simulation_file] = obj.get_path(obj.config);
            simulation_path = [simulation_dir,simulation_file];
            if ~exist(simulation_dir,'dir')
                mkdir(simulation_dir)
            end
            simulation.config.drogued_status = obj.config.drogued_status;
            simulation.config.dx = obj.config.dx;
            simulation.config.dt = obj.config.dt;
            simulation.grid = obj.grid;
            simulation.time = obj.time;
            simulation.tracer = obj.tracer;
            simulation.config.days_in_year = obj.config.days_in_year;
            save(simulation_path,'simulation','-v7.3');
        end
    end
       
    methods(Static)
        function simulation = load(config)
            % Loads existing simulation from a .mat file (static method).
            %
            % Input arguments:
            % - drogued_status: string that can be: "all", "drogued" or "undrogued"
            % - dx [degrees]: grid size
            % - dt [days]: time step
            % - sinks[logical]: 0: sinks are removed, 1: sinks are fixed
            % - run_time [years]: simulation run time
            % - output_interval [integer]: interval for which simulation
            %   output is stored
            %
            % Output:
            % - simulation: data structure containing simulation result
            %
            [simulation_dir,simulation_file] = TmModel.get_path(config);
            if exist([simulation_dir,simulation_file],'file')
                load([simulation_dir,simulation_file]);
            else
                error('TmModel:fileNotFound','Simulation does not yet exist. Please run it first');
            end
        end
        
        function [simulation_dir,simulation_file] = get_path(config)
            % Gets simulation directory and simulation filename (static method).
            %
            % Input arguments:
            % - drogued_status: string that can be: "all", "drogued" or "undrogued"
            % - initial_type: string that describes initial condition type,
            %   existing options are: "global_uniform", "uniform",
            %   "jambeck2015", "lebreton2017"
            % - dx [degrees]: grid size
            % - dt [days]: time step
            % - sinks [logical]: 0: sinks are removed, 1: sinks are fixed
            % - run_time [years]: simulation run time
            % - output_interval [integer]: interval for which simulation
            %   output is stored
            %
            % Output:
            % - simulation_dir: directory where .mat file with simulation
            %   is stored
            % - simulation_file: filename where .mat file with simulation
            %   is stored
            %
            [tm_dir,~] = TransportMatrix.get_path(config);
            if strcmpi(config.initial.type,'global_uniform')
                initial_description = '_gu';
            elseif strcmpi(initial_type,'uniform')
                initial_description = '_u';
            elseif strcmpi(initial_type,'jambeck2015')
                initial_description = '_jb2015';
            elseif strcmpi(initial_type,'lebreton2017')
                initial_description = '_lb2017';
            else
                error('Unknown initial_type.');
            end
            run_time_description = ['rt',num2str(config.run_time)];
            output_interval_description = ['_oi',num2str(config.output_interval)];            
            simulation_dir = [tm_dir,run_time_description,output_interval_description,initial_description,'/'];
            simulation_file = 'sim.mat';
        end
        
        function [time_index,passed_time] = get_time_index(simulation,requested_time)
            % Gets index of requested_time in simulation time (static method).
            % If the requested time is not available in the simulation, the
            % nearest time to the requested time is returned.
            %
            % Input arguments:
            % - simulation: data structure containing simulation result.
            % - requested_time [years]: time for which index in simulation
            %   time is returned, unit is years after simulation start
            %
            % Output:
            % - time_index: index of simulation time that matches or is
            %   closest to the requested time
            % - passed_time: time that has passed since the start of the
            %   simulation (this is the same as requested_time, unless the
            %   requested_time is not available and the closest time is
            %   returned)
            %
            time_since_start = simulation.time-simulation.time(1);
            dt = time_since_start-requested_time*simulation.config.days_in_year;
            closest_time = find(abs(dt)==min(abs(dt)));
            time_index = closest_time(1);
            passed_time = time_since_start(time_index)/simulation.config.days_in_year;
        end
        
        function tracer_2d = convert_1d_to_2d(tracer_1d,grid)
            % Converts 1d tracer [loc] to 2d [lat,lon] (static method). If
            % tracer has a time component as well, use convert_2d_to_3d
            % instead.
            %
            % Input arguments:
            % - tracer_1d [loc]: vector containing tracer concentration
            % - grid: GlobalGrid object
            %
            % Output:
            % - tracer_2d [lat,lon]: 2d matrix containing tracer
            %   concentration, dimensions match grid.lon and grid.lat            
            %            
            tracer_2d(:,:) = vec2mat(tracer_1d,grid.lon_size);
        end
        
        function tracer_3d = convert_2d_to_3d(tracer_2d,grid)
            % Converts 2d simulation result [time,loc] to 3d [lat,lon,time]
            % (static method).
            %
            % Input arguments:
            % - tracer_2d [time,loc]: matrix containing tracer concentration
            % - grid: GlobalGrid object
            %
            % Output:
            % - tracer_3d [lat,lon,time]: 3d matrix containing tracer
            %   concentration, dimensions match grid.lon and grid.lat
            %   and original tracer_2d time dimensions
            %
            temporal_size = size(tracer_2d,1);
            tracer_3d = nan(grid.lat_size,grid.lon_size,temporal_size);
            for i = 1:temporal_size
                tracer_3d(:,:,i) = vec2mat(tracer_2d(i,:),grid.lon_size);
            end
        end
    end
end