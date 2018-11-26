classdef GlobalGrid < handle
    properties
        resolution
        lon_min
        lon_max
        lat_min
        lat_max
        lon_size
        lat_size
        lon
        lat
    end
    
    methods
        function obj = GlobalGrid(resolution)
            % Creates a global grid with square grid cells.
            %
            % Input arguments:
            % - resolution: (degrees) resolution used to build grid cells.
            %
            % Output:
            % - GlobalGrid object
            %
            obj.resolution = resolution;
            range_lon = [-180,180];
            range_lat = [-90,90];
            obj.lon_max = max(range_lon);
            obj.lon_min = min(range_lon);
            obj.lat_max = max(range_lat);
            obj.lat_min = min(range_lat);
            obj.lon_size = abs(obj.lon_max-obj.lon_min)*1/obj.resolution;
            obj.lat_size = abs(obj.lat_max-obj.lat_min)*1/obj.resolution;
            obj.lon = [obj.lon_min:obj.resolution:obj.lon_max-obj.resolution];
            obj.lat = [obj.lat_min:obj.resolution:obj.lat_max-obj.resolution];            
        end
        
        function [lon_index,lat_index] = get_index(obj,lon,lat)
            % Gets index of a longitude, latitude location in a grid.
            %
            % Input arguments:
            % - lon: longitude of a location
            % - lat: latitude of a location
            %
            % Output:
            % - lon_index: x-index in grid of longitude location
            % - lat_index: y-index in grid of latitude location
            %
            lon_index = floor(lon*1/obj.resolution)-obj.lon_min*1/obj.resolution+1;
            l_index_lon_over = lon_index==abs(obj.lon_max-obj.lon_min)*1/obj.resolution+1;
            if any(l_index_lon_over)
               lon_index(l_index_lon_over) = 1;
            end
            lat_index = floor(lat*1/obj.resolution)-obj.lat_min*1/obj.resolution+1;
            l_index_lat_over = lat_index==abs(obj.lat_max-obj.lat_min)*1/obj.resolution;
            if any(l_index_lat_over)
                lat_index(l_index_lat_over) = 1;
            end
        end
       
    end
end
