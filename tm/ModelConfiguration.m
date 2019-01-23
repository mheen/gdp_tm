classdef ModelConfiguration < handle
    properties
        description
        drogued_status
        dx
        dt        
        run_time
        output_interval
        sinks
        normalise
        initial
        days_in_year
        specify_time
    end
    
    methods        
        function obj = ModelConfiguration(varargin)
            % Gets settings for the transport matrix model simulation.
            % Default settings are loaded first and then overwritten by
            % settings contained in an *.ini file. The set-up of the
            % *.ini file is explained in the README file. An example
            % "config.ini" file with instructions is also included.
            %
            % Input arguments:
            % - input_path: path to the *.ini file that contains model
            %   settings
            %
            % Requirements:
            % This function uses "IniConfig.m" distributed through
            % MathWorks File Exchange by Evgeny Prilepin aka Iroln. This
            % script and the accompanying license has been copied into the
            % "../external_functions/" folder for user convenience.
            %
            obj.get_defaults()
            obj.read_config(varargin)            
        end
        
        function read_config(obj,varargin)
            % Reads model settings from either an *.ini file or from
            % parameters given as function input (variable number of
            % inputs).
            %
            if length(varargin{:}) == 1
                input_path = char(varargin{1});
                if strcmpi(input_path(end-3:end),'.ini')
                    obj.read_config_from_ini(input_path)
                else
                    error('Configuration input file must be a .ini file.');
                end
            else
                obj.read_config_parameters(varargin{:});
            end
        end
        
        function read_config_from_ini(obj,input_path)
            % Read model settings from the *.ini file.
            %
            % Input arguments:
            % - input_path: path to the *.ini file
            %
            config = IniConfig();
            config.ReadFile(input_path);
            obj.description = config.GetSections{1};
            config_keys = config.GetKeys(obj.description);
            
            for i = 1:length(config_keys)
                key_part = strfind(config_keys{i},'.');
                if ~isempty(key_part)
                    first_key_name = config_keys{i}(1:key_part-1);
                    second_key_name = config_keys{i}(key_part+1:end);
                    obj.(first_key_name).(second_key_name) = config.GetValues(obj.description,config_keys{i});
                else
                    obj.(config_keys{i}) = config.GetValues(obj.description,config_keys{i});
                end
            end            
        end
        
        function read_config_parameters(obj,varargin)
            % Reads model settings from input parameters.
            %
            for n = 1:2:length(varargin{:})-1
                if strcmpi(varargin{1}{n},'drogued_status')
                    obj.drogued_status = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'dx')
                    obj.dx = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'dt')
                    obj.dt = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'run_time')
                    obj.run_time = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'output_interval')
                    obj.output_interval = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'sinks')
                    obj.sinks = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'normalise')
                    obj.normalise = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'initial.type')
                    obj.initial.type = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'initial.lon')
                    obj.initial.lon = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'initial.lat')
                    obj.initial.lat = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'initial.tracer_per_grid_cell')
                    obj.initial.tracer_per_grid_cell = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'initial.years_to_add_sources')
                    obj.initial.years_to_add_sources = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'days_in_year')
                    obj.days_in_year = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'specify_time.type')
                    obj.specify_time.type = varargin{1}{n+1};
                elseif strcmpi(varargin{1}{n},'specify_time.times')
                    obj.specify_time.times = varargin{1}{n+1};
                end
            end            
        end
        
        function get_defaults(obj)
            % Loads default model settings.
            %
            obj.description = '[tm model configuration]';
            obj.drogued_status = 'all'; % 'all', 'drogued' or 'undrogued'
            obj.dx = 1; % degrees
            obj.dt = 60; % days            
            obj.run_time = 50; % years
            obj.output_interval = 1;
            obj.sinks = 0;
            obj.normalise = 0;
            obj.initial.type = 'global_uniform';
            obj.initial.lon = [];
            obj.initial.lat = [];
            obj.initial.tracer_per_grid_cell = 1;
            obj.initial.years_to_add_sources = [];            
            obj.days_in_year = 360; % this needs to be a multiple of dt
            obj.specify_time.type = 'all';
            obj.specify_time.times = [];
        end
        
        function check_normalisation(obj)
            % Checks if normalisation is set correctly when using an
            % initial condition that injects sources. Normalisation is set
            % to 1 if this is the case. Note that this overrules the
            % normalisation settings in the config.ini file.
            %
            if strcmpi(obj.initial.type,'jambeck2015') || strcmpi(obj.initial.type,'lebreton2017')
                obj.normalise = 1;
            end
        end
        
        function write_config(obj,output_dir)
            % Writes model configuration settings to a config.ini file.
            % This is useful when running several different simulations,
            % the corresponding config.ini file can then be written to the
            % simulation output directory.
            %
            % Input arguments:
            % - output_dir: directory where config.ini file is written to
            %
            output_path = [output_dir,'config.ini'];
            config = IniConfig();
            config.AddSections(obj.description);
            config.AddKeys(obj.description,'drogued_status',obj.drogued_status);
            config.AddKeys(obj.description,'dx',obj.dx);
            config.AddKeys(obj.description,'dt',obj.dt);           
            config.AddKeys(obj.description,'run_time',obj.run_time);
            config.AddKeys(obj.description,'output_interval',obj.output_interval);
            config.AddKeys(obj.description,'sinks',obj.sinks);
            config.AddKeys(obj.description,'normalise',obj.normalise);
            config.AddKeys(obj.description,'initial.type',obj.initial.type);
            config.AddKeys(obj.description,'initial.lon',obj.initial.lon);
            config.AddKeys(obj.description,'initial.lat',obj.initial.lat);
            config.AddKeys(obj.description,'initial.tracer_per_grid_cell',obj.initial.tracer_per_grid_cell);
            config.AddKeys(obj.description,'initial.years_to_add_sources',obj.initial.years_to_add_sources);
            config.AddKeys(obj.description,'days_in_year',obj.days_in_year);
            config.AddKeys(obj.description,'specify_time.type',obj.specify_time.type);
            config.AddKeys(obj.description,'specify_time.times',obj.specify_time.times);
            config.WriteFile(output_path);          
        end                
    end    
end