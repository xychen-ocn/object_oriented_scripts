% ======================================================================= %
% class: ATOMIC_GOES16
% purpose: process ATOMIC data:
% Author: Xuanyu Chen 
% Date: v0 @ 09/29/2021
% ======================================================================= %

classdef ATOMIC_LagrangianTraj
    properties
        data struct
        TrajModel
        SimInputDir = '/Users/xchen/Documents/MATLAB/shallow_convection/generate_lagtraj'
        SimOutputDir 
        SimID
        SimDuration
        TrajLevs
        ATOMIC_platform
   
    end
    
    methods
        %%%% ===================   function #0     =================== %%%%
        % Purpose: initialize the object:
        function obj = ATOMIC__LagrangianTraj(inval)
            if nargin > 0 
                if isstruct(inval)
                    fieldns = fieldnames(inval);
                    for i = 1:length(fieldns)
                        FN = fieldns{i};
                        obj.(FN) = inval.(FN);
                    end
              
                end
            else
                obj.data = struct();  % initialize with an empty structure.
                obj.TrajModel = 'LagTraj';
                obj.SimOutputDir = '';
                obj.SimID = '';
                obj.SimDuration = 6;
                obj.TrajLevs = 700;
                obj.ATOMIC_platform = '';
            end
            
        end
        
        
        %%%% ===================   function #1     =================== %%%%
        % Purpose: generate input for HYSPLIT model to generate a
        % trajectory
        function generate_input_for_trajectory_simulation(obj)
            % The input I need will be:
            % 1. the location of the start point
            % 2. the start time of the simulation
            % 3. the height of the trajectory
            % 4. the duration of simulation (equally upstream and down
            % stream)
            % This information will be written into a text file as input
            % for the HYSPLIT model?
            %
            % I can also write this into a yaml file if I am going to use
            % the lagtraj toolbox.
            
            dataIn = obj.data;
            for iseg = 1:length(dataIn)
                origin.lon(iseg) = mean(dataIn(iseg).lon);
                origin.lat(iseg) = mean(dataIn(iseg).lat);
                simparam.T0_UTC(iseg) = mean(dataIn(iseg).time);      % UTC
                simparam.duration = obj.SimDuration;                     % same for both backward and forward traj
                simparam.height = obj.TrajLevs;
            end
            
            % make input directory to save these input files:
            inp_svdir = [obj.SimInputDir filesep obj.TrajModel filesep 'inp_yaml_files'];
            if ~exist(inp_svdir)
                mkdir(inp_svdir)
            end
            
            % write output in corresponding format
            if strcmpi(obj.TrajModel, 'HYSPLIT')
               
                write_HYSPLIT_input(origin, simparam)
                
            elseif strcmpi(obj.TrajModel, 'lagtraj')
                
                write_yaml_input(origin, simparam)
                
            end
            
            
            %%%%%%%%%%%%%%%% Nested function 01
            function write_HYSPLIT_input()
                
                % write file into the input directory 
                
            end
            
            
            %%%%%%%%%%%%%%%% Nested function 02
            function write_yaml_input(origin, simparam)
                % this format is much easier or simpler.
                
                % 1. figure out how many files we are writing out: = num of
                % simulation origins.
                num_yaml_files = length(origin.lon);
                fmt1 = 'yyyy-mm-ddTHH:MM';
                
                % 2. write out these files in batches:
                for i = 1:num_yaml_files
                    dateIDstr = datestr(simparam.T0_UTC(i), 'yyyymmdd_HHMM');
                    filename = ['ATOMIC_' obj.ATOMIC_platform '_' ...
                                dateIDstr '_seg' num2str(i,'%2.2i')...
                                '.yaml'];
                    absFN = [inp_svdir filesep filename];
                    
                    fid = fopen(absFN,'w');
                    
                    % construct parameter fields in the yaml file:
                    
                    fprintf(fid, '%s\n', 'trajectory_type : lagrangian');
                    fprintf(fid, '%s\n', 'version         : 1.0.0');
                    fprintf(fid, '%s\n', 'velocity_method : single_height_level');
                    fprintf(fid, '%s %5.1f\n', 'velocity_method_height :', simparam.height);
                    fprintf(fid, '%s\n', 'domain          : eurec4a_northwest_atlantic');
                    fprintf(fid, '%s %7.3f\n', 'lat_origin      :', origin.lat(i));
                    fprintf(fid, '%s %7.3f\n', 'lon_origin      :', origin.lon(i));
                    fprintf(fid, '%s\n', ['datetime_origin : "', ...
                          datestr(simparam.T0_UTC(i),fmt1), '"']);
                    fprintf(fid, '%s%iH\n', 'backward_duration : PT', simparam.duration);
                    fprintf(fid, '%s%iH\n', 'forward_duration : PT', simparam.duration);
                    fprintf(fid, '%s\n', 'timestep        : PT0.25H');   % 15min time step.
                    
                    fclose(fid);
                    disp(['-> Finished writing input yaml file:' filename]);
                end
            end
            
                
            
            
        end
        
        %%%% ===================   function #2     =================== %%%%
        % Purpose: get output from the lagtraj toolbox
        function lagtraj = get_output_from_trajectory_simulation(obj)
            % I need to know the output format first.. (could be ASCII)     
            if strcmp(obj.TrajModel, 'LagTraj')
                % read everything from a netCDF file:
                outdir = [obj.SimOutputDir filesep obj.TrajModel filesep 'out_netCDF'];
                files = dir([outdir filesep 'ATOMIC_' obj.ATOMIC_platform '*.nc']);
                filenames = {files.name};
                NumFiles = length(filenames);
                
                %lagtraj = struct();
                for iseg= 1:NumFiles
                    disp(['--> reading ' filenames{iseg}]);
                    absFN = [outdir filesep filenames{iseg}];
                    lagtraj(iseg) = read_netCDF_into_matlab_structure(absFN);
                end
                
                % each structure will contain the following fields:
                % . time
                % . lat
                % . lon
                % . u_traj
                % . v_traj
                % . origin_lat
                % . origin_lon
                % . origin_datetime
                
            elseif strcmp(obj.TrajModel, 'HYSPLIT')
                lagtraj = get_output_from_HYSPLIT(obj, opt_args);
            end
             
        end
        
        
        
        
        %%%% ===================   function #3     =================== %%%%
        % Purpose: get output from HYSPLIT model 
        function lagtraj = get_output_from_HYSPLIT(obj, fileIDs)
            % I need to know the output format first.. (could be ASCII)
            
                % collect the output from noaa.gov.jhypubout..
                %urlbase = 'https://www.ready.noaa.gov/hypubout';
                urlbase = obj.SimOutputDir;
                
                lagtraj = struct();
                for iseg = 1:length(fileIDs)
                    btrajID = fileIDs(iseg).btraj;
                    ftrajID = fileIDs(iseg).ftraj;
                    file_backward = [urlbase filesep 'tdump.' btrajID '.txt'];
                    file_forward = [urlbase filesep 'tdump.' ftrajID '.txt'];
                    
                    % read these two files:
                    bw_endpoints = read_HYSPLIT_weboutput(file_backward);
                    fw_endpoints = read_HYSPLIT_weboutput(file_forward);
                    
                    
                    
                    % combine backward traj and forward traj:
                    fields = fieldnames(bw_endpoints);
                    for i = 1:length(fields)
                        FN = fields{i};
                        all_endpnts.(FN) = unique([bw_endpoints.(FN); fw_endpoints.(FN)]);
                    end
                    
                    
                    lagtraj(iseg).lon = all_endpnts.lon;
                    lagtraj(iseg).lat = all_endpnts.lat;
                    lagtraj(iseg).time = all_endpnts.time;
                end
            
            %%%%%%%%%%%%%%%% Nested function
            function read_HYSPLIT_weboutput(fileurl)
                data = webread(fileurl);
                % check how Ali read this type of data and break it down.
                
                % put data into categories according to the input object.
                % that will help to determine how many height levels are
                % used,
                num_levs = length(obj.TrajLevs);
                
                
                
            end
            
            
        end
        

        

        

    end
    
end