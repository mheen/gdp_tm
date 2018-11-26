classdef AccumulationRegions < handle
    properties
        tracer_threshold % minimum tracer amount per grid cell to define an accumulation
        area_threshold % minimum number of adjacent grid cells that make up an accumulation
        time % simulation time
        grid % GlobalGrid object containing information on transport matrix grid
        global_tracer % simulated tracer concentration
        names % names of accumulation regions
        outlines % outlines of accumulation regions
        tracer % tracer contained in accumulation regions
        total_tracer % total tracer in simulation
        max_tracer % maximum tracer values and locations within accumulation regions
    end
    
    methods
        function obj = AccumulationRegions(simulation,config,tracer_threshold,area_threshold)
            % Accumulation regions determined from simulation output.
            % Defined as regions with at least area_threshold adjacent grid
            % cells that contain a tracer concentration above the
            % tracer_threshold value. Outlines of all regions that satisfy
            % these condition are found first. Then, specific areas in
            % regions of interest (subtropical accumulations) are found.
            % The outlines for these regions are smoothed and the total
            % number of tracer in the accumulation regions is determined.
            %
            % Input arguments:
            % - simulation: data structure containing simulation result, created by
            %   running a transport matrix simulation using "TmModel"
            % - tracer_threshold: tracer amount per grid cell that is set as a minimum
            %   value to determine accumulation regions
            % - area_threshold: minimum number of adjacent grid cells that make up an
            %   accumulation region
            %
            obj.tracer_threshold = tracer_threshold;
            obj.area_threshold = area_threshold;
            obj.time = simulation.time;
            obj.grid = simulation.grid;
            obj.global_tracer = simulation.tracer;
            
            obj.get_all_outlines()
            obj.get_specific_outlines()
            obj.get_tracer_in_accumulations()
            obj.save(config)
        end
        
        function get_tracer_in_accumulations(obj)
            % Gets the total amount of tracer in specific accumulation
            % regions as well as the maximum value and location.
            %
            for t = 1:length(obj.time)
                for i = 1:length(obj.names)
                    l_region = ~isnan(obj.outlines.(obj.names{i}){t}.matrix);
                    global_tracer_t = obj.global_tracer(:,:,t);
                    obj.tracer.(obj.names{i}){t} = sum(sum(global_tracer_t(l_region)));
                    % get maximum tracer value and location
                    max_value = max(max(global_tracer_t(l_region)));
                    if ~isempty(max_value)
                        obj.max_tracer.(obj.names{i}).value{t} = max_value;
                        obj.max_tracer.(obj.names{i}).index{t} = find(global_tracer_t==max_value);
                        [i_max,j_max] = ind2sub(size(global_tracer_t),obj.max_tracer.(obj.names{i}).index{t});
                        obj.max_tracer.(obj.names{i}).i{t} = i_max;
                        obj.max_tracer.(obj.names{i}).j{t} = j_max;
                        [lon_max,lat_max] = obj.grid.get_index(i_max,j_max);
                        obj.max_tracer.(obj.names{i}).lon{t} = lon_max;
                        obj.max_tracer.(obj.names{i}).lat{t} = lat_max;
                    end
                end
            end
            obj.total_tracer = sum(sum(obj.global_tracer(:,:,1)));
        end
        
        function get_specific_outlines(obj)
            % Extracts outlines of specific regions from all outlines.
            %
            region_boxes.np = obj.get_accumulation_box('np',obj.grid);
            region_boxes.sp = obj.get_accumulation_box('sp',obj.grid);
            region_boxes.na = obj.get_accumulation_box('na',obj.grid);
            region_boxes.sa = obj.get_accumulation_box('sa',obj.grid);
            region_boxes.si = obj.get_accumulation_box('si',obj.grid);
            obj.names = fieldnames(region_boxes);
            
            for t = 1:length(obj.time)
                for i = 1:length(obj.names)
                    obj.outlines.(obj.names{i}){t}.matrix = nan(size(obj.outlines.all{t}.matrix));
                    i_lon = [region_boxes.(obj.names{i}).i_lon(1):region_boxes.(obj.names{i}).i_lon(2)];
                    i_lat = [region_boxes.(obj.names{i}).i_lat(1):region_boxes.(obj.names{i}).i_lat(2)];
                    n_count = 1;
                    if isfield(obj.outlines.all{t},'labels')
                        for j = 1:length(obj.outlines.all{t}.labels)
                            i_area = obj.outlines.all{t}.matrix==obj.outlines.all{t}.labels(j);
                            region = i_area(i_lat,i_lon);
                            if sum(sum(region)) > obj.area_threshold
                                x = obj.outlines.all{t}.i{j};
                                y = obj.outlines.all{t}.j{j};
                                [smooth_x,smooth_y] = obj.smooth_boundaries(x,y);
                                mask_matrix = poly2mask(smooth_x,smooth_y,size(obj.outlines.all{t}.matrix,1),size(obj.outlines.all{t}.matrix,2));
                                smooth_matrix = obj.outlines.all{t}.matrix.*mask_matrix;
                                smooth_matrix(smooth_matrix==0) = NaN;
                                if any(any(smooth_matrix(i_area)))
                                    obj.outlines.(obj.names{i}){t}.i{n_count} = smooth_x;
                                    obj.outlines.(obj.names{i}){t}.j{n_count} = smooth_y;
                                    obj.outlines.(obj.names{i}){t}.lon{n_count} = obj.grid.lon(smooth_x);
                                    obj.outlines.(obj.names{i}){t}.lat{n_count} = obj.grid.lat(smooth_y);
                                    obj.outlines.(obj.names{i}){t}.matrix(i_area) = smooth_matrix(i_area);
                                    obj.outlines.(obj.names{i}){t}.labels(n_count) = obj.outlines.all{t}.labels(j);
                                    n_count = n_count+1;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function get_all_outlines(obj)
            % Determines outlines of all accumulation regions that meet the
            % requirements of:
            % 1. The tracer concentration is above the tracer_threshold.
            % 2. The area where tracer concentration is above the
            %    tracer_threshold is larger than area_threshold grid cells.
            %
            for t = 1:length(obj.time)
                l_tracer_above_threshold = obj.global_tracer(:,:,t) > obj.tracer_threshold;
                % l_tracer_above_threshold is used as a black/white (binary)
                % image, around which boundaries are drawn
                [boundaries,label_matrix] = bwboundaries(l_tracer_above_threshold,'noholes');
                boundary_labels = unique(label_matrix);
                boundary_labels = boundary_labels(2:end); % start at 2 because boundary_labels contains 0 label
                area_in_boundaries = zeros(length(boundary_labels),1);
                for b = 1:length(boundary_labels)
                    area_in_boundaries(b) = sum(sum(label_matrix==boundary_labels(b)));
                end
                % get outlines with area > area_threshold
                n_count = 1;
                obj.outlines.all{t}.matrix = nan(size(label_matrix));
                for b = 1:length(boundaries)
                    if area_in_boundaries(b) > obj.area_threshold
                        l_label = label_matrix == boundary_labels(b);
                        obj.outlines.all{t}.matrix(l_label) = label_matrix(l_label);
                        obj.outlines.all{t}.labels(n_count) = boundary_labels(b);
                        obj.outlines.all{t}.i{n_count} = boundaries{b}(:,2);
                        obj.outlines.all{t}.j{n_count} = boundaries{b}(:,1);
                        obj.outlines.all{t}.lon{n_count} = obj.grid.lon(boundaries{b}(:,2));
                        obj.outlines.all{t}.lat{n_count} = obj.grid.lat(boundaries{b}(:,1));
                        n_count = n_count+1;
                    end
                end
            end
        end
        
        function save(obj,config)
            % Writes the accumulation regions to a data structure and saves
            % to a .mat file. The .mat file is saved to the directory that
            % also contains the relevant simulation.
            %
            [ac_dir,ac_file] = obj.get_path(config,obj.tracer_threshold,obj.area_threshold);
            ac_path = [ac_dir,ac_file];
            if ~exist(ac_dir,'dir')
                mkdir(ac_dir)
            end
            accumulation.tracer_threshold = obj.tracer_threshold;
            accumulation.area_threshold = obj.area_threshold;
            accumulation.time = obj.time;
            accumulation.grid = obj.grid;
            accumulation.names = obj.names;
            accumulation.outlines = obj.outlines;
            accumulation.tracer = obj.tracer;
            accumulation.total_tracer = obj.total_tracer;
            accumulation.max_tracer = obj.max_tracer;
            save(ac_path,'accumulation','-v7.3');
        end        
    end
    
    
    methods(Static)
        function [smooth_x,smooth_y] = smooth_boundaries(x,y)
            % Smoothes outlines of accumulation regions by getting rid of
            % crossing minor polygons.
            %
            polygon = [x,y];
            [~,~,ic] = unique(polygon,'rows','stable');
            [n,bin] = histc(ic,unique(ic));
            multiple = find(n>1);
            for i = 1:length(multiple)
                index_duplicates = find(bin==multiple(i));
                internal_polygon = polygon(index_duplicates(1)+1:index_duplicates(2)-1,:);
                [~,~,internal_ic] = unique(internal_polygon,'rows','stable');
                [internal_n,internal_bin] = histc(internal_ic,unique(internal_ic));
                internal_multiple = find(internal_n>1);
                for j = 1:length(internal_multiple)
                    internal_index_duplicates = find(internal_bin==internal_multiple(j));
                    start_remove = internal_index_duplicates(1)+1;
                    end_remove = internal_index_duplicates(end);
                    internal_polygon(start_remove:end_remove,:) = NaN;
                end
                polygon_number = ['n',num2str(i)];
                internal_polygon_no_nan = internal_polygon(~isnan(internal_polygon(:,1)),:);
                start_point = polygon(index_duplicates(1),:);
                end_point = polygon(index_duplicates(2),:);
                smooth_polygons.(polygon_number) = [start_point;internal_polygon_no_nan;end_point];
                polygon_area(i) = polyarea(smooth_polygons.(polygon_number)(:,1),smooth_polygons.(polygon_number)(:,2));
                clear internal_polygon
            end
            [~,n_smooth_polygon] = max(polygon_area);
            polygon_number = ['n',num2str(n_smooth_polygon)];
            smooth_polygon = smooth_polygons.(polygon_number);
            smooth_x = smooth_polygon(:,1);
            smooth_y = smooth_polygon(:,2);
        end
        
        function box = get_accumulation_box(region,grid)
            % Returns a box with a longitude and latitude range where to
            % look for a specific accumulation region, as well as the
            % corresponding indices in the grid.
            %
            % Input arguments:
            % - region: string with name of specific accumulation region
            %   valid options are: 'np','sp','na','sa','si' for the
            %   subtropical accumulation regions, and 'barents',
            %   'baybengal', 'eastaus' for other regions.
            %
            % Output:
            % - box: datastructure containing
            %   * lon: array with minimum and maximum longitude range
            %   * lat: array with minimum and maximum latitude range
            %   * i_lon: array with indices in grid corresponding to lon
            %   * i_lat: array with indices in grid corresponding to lat
            %
            if strcmpi(region,'np')
                box.lon = [-180 -120];
                box.lat = [22.5 37.5];
            elseif strcmpi(region,'na')
                box.lon = [-60 -40];
                box.lat = [22.5 37.5];
            elseif strcmpi(region,'sp')
                box.lon = [-180 -85];
                box.lat = [-37.5 -22.5];
            elseif strcmpi(region,'sa')
                box.lon = [-30 0];
                box.lat = [-37.5 -22.5];
            elseif strcmpi(region,'si')
                box.lon = [50 80];
                box.lat = [-37.5 -22.5];
            elseif strcmpi(region,'barents')
                box.lon = [0 60];
                box.lat = [60 90];
            elseif strcmpi(region,'baybengal')
                box.lon = [60 120];
                box.lat = [0 30];
            elseif strcmpi(region,'eastaus')
                box.lon = [120 180];
                box.lat = [-40 -20];
            else
                error(['Box requested for unknown accumulation region. Valid options:',...
                    'np,sp,na,sa,si,barents,baybengal,eastaus'])
            end
            [box.i_lon(1),box.i_lat(1)] = grid.get_index(box.lon(1),box.lat(1));
            [box.i_lon(2),box.i_lat(2)] = grid.get_index(box.lon(2),box.lat(2));
        end
        
        function [ac_dir,ac_file] = get_path(config,tracer_threshold,area_threshold)
            % Gets directory and filename to save accumulation regions to
            % (saved to the same folder as the relevant simulation).
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
            % - ac_dir: directory where .mat file with accumulation regions
            %   is stored
            % - ac_file: filename where .mat file with accumulation regions
            %   is stored
            %
            [ac_dir,~] = TmModel.get_path(config);
            tracer_description = ['_tt',num2str(tracer_threshold)];
            area_description = ['_at',num2str(area_threshold)];
            ac_file = ['ac',tracer_description,area_description,'.mat'];
        end
        
        function accumulation = load(config,tracer_threshold,area_threshold)
            [ac_dir,ac_file] = AccumulationRegions.get_path(config,tracer_threshold,area_threshold);            
            if exist([ac_dir,ac_file])
                load([ac_dir,ac_file])
            else
                simulation = TmModel.load(config);
                accumulation = AccumulationRegions(simulation,config,tracer_threshold,area_threshold);
            end
        end
        
    end
end
