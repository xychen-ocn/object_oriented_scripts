% ======================================================================= %
% class: ATOMIC_GOES16
% purpose: process ATOMIC data:
% Author: Xuanyu Chen 
% Date: v0 @ 09/29/2021
% ======================================================================= %

classdef ATOMIC_GOES16
    properties
        data struct
        SateSource
        SateDir
        SateChannel
        type 
        DateLST
        LocSaveDir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/saved_SateData';
        access = ''
   
    end
    
    methods
        
        %%%% ===================   function #0     =================== %%%%
        % Purpose: initialize the object:
        function obj = ATOMIC_GOES16(inval)
            if nargin > 0 
                if isstruct(inval)
                    fieldns = fieldnames(inval);
                    for i = 1:length(fieldns)
                        FN = fieldns{i};
                        obj.(FN) = inval.(FN);
                    end
                    
                    % for fields that are not included in the inval,
                    % initialize with empty :
              
                end
            else
                obj.data = struct();  % initialize with an empty structure.
                obj.SateSource = '';
                obj.SateDir = '';
                obj.SateChannel = '';
                obj.type = '';
                obj.DateLST = '';
                obj.access = '';
            end
            
        end
        
        
        
        %%%% ===================   function #1     =================== %%%%
        % Purpose: get data from GOES16 based on input
        function data_obj = extract_GOES16_data(obj, SateChannel, box_width, time_window)
            % 0. initailize the data obj:
            data_obj = ATOMIC_GOES16();
            data_obj.SateSource = obj.SateSource;
            data_obj.SateDir = obj.SateDir;
            data_obj.SateChannel = SateChannel;
            data_obj.type = 'SateData';
            data_obj.DateLST = obj.DateLST;
            
            % 1. determine source (product name, IR or VIS)
            % .  each source will correspond to slightly different file
            % names and subfunction to use to extract data effectively.
            % filename will need to be constructed dynamcially at each
            % time instance.
            
            % determine channel:
            if strcmpi(SateChannel, 'VIS')
                channelID='02';
                
            else
                channelID='13';
                
            end
            
            % 2. find out number of segment in the obj.data
            num_segs = length(obj.data);
            for iseg = 1:num_segs
                seg_time = obj.data(iseg).time;
                trange =  mean(seg_time)-0.5*time_window/24:tres: mean(seg_time)+0.5*time_window/24;
                Nt = length(trange);
                
                for it =1:Nt
                    tnow = trange(it);
                    disp(['** tnow:' datestr(tnow)]);
                    % 1. construct filename dynamically:
                    absFN = construct_filename(channelID, tnow);
                    
                    % 2. read and extract subset of data:  a region that is large enough to cover the
                    % lagrangian trajectory and save the data
                    datasub_tmp = read_subset_of_data(absFN, scene_center, tnow, box_width);
                    
                    datasub(iseg).lon(:,it) = squeeze(datasub_tmp.lon);
                    datasub(iseg).lat(:,it) = squeeze(datasub_tmp.lat);
                    datasub(iseg).time(:,it) = tnow;
                    datasub(iseg).values(:,:,it) = squeeze(datasub_tmp.values);
                    
                    
                end
                
            end
            
            % 3. save data:
            matFN = [obj.DateLST 'LST_centered_on_' obj.type '_boxwidth' num2str(box_width,'%2.2i') '.mat'];
            
            % establish a subfolder if it doesn't exist.
            svdir = [obj.LocSaveDir filesep obj.SateSource];
            if ~exist(svdir,'dir')
                mkdir(svdir);
            end
            
            data_obj.data = datasub;
            data_obj.access = [svdir filesep matFN];
            
            save([svdir filesep matFN], 'data_obj');
            
            
            % datasub contains all segments:
            
            
            
            
            
            %%%%%%%%%%%%%%%%%% Nested function 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function absFN = construct_filename(channelID, tnow)
                format1 = 'yyyy_mm_dd';
                format2 = 'yyyymmdd';
                
                subfolderN = datestr(floor(tnow),format1);
                path2data = [obj.SateDir filesep obj.SateSource filesep '2020' filesep subfolderN];
                
                doy = floor(tnow) - datenum(2020,1,1)+1; % day of year;
                HH = datestr(tnow,'HH');
                MM = datestr(tnow,'MM');
                
                MMr = round(str2double(MM)/10)*10;
                if MMr ==60
                    HHstr = num2str(str2double(HH)+1, '%2.2i');
                    MMr = 0;
                else
                    HHstr = HH;
                end
                MMr_str = num2str(MMr,'%2.2i');
                
                if strcmp(obj.SateSource, '0.5km_01min')
                    NCFN = ['GOES_' channelID '_8N-18N-62W-50W_' datestr(floor(tnow),format2) '.nc'];
                    absFN = [path2data filesep NCFN];
                    
                elseif strcmp(obj.SateSource, '2km_10min')
                    
                    % read data off the web catalog:
                    url0 = [path2data filesep 'catalog.html'];
                    url_thredds  = replace(url0, 'dodsC','catalog');
                    
                    url_backup = ['https://observations.ipsl.fr/aeris/eurec4a-data/SATELLITES/GOES-E/2km_10min/2020/' subfolderN filesep];  % backup:
                    [tmp, url_stat] = urlread(url0);
                    
                    if url_stat == 0
                        % invalid url, use the new url instead:
                        url_valid = url_backup;
                    else
                        url_valid = url0;
                    end
                    
                    fileid = ['2020' num2str(doy,'%3.3i') HHstr MMr_str]; % YYYYJJJHHMMSSZ  (4digit year, 3digit day of year, hour, minute, second, tenth second
                    disp(['**dataID:' fileid]);
                    
                    fileid_full = get_fileID_from_thredds_catalog(url_valid, fileid);
                    
                    NCFN = ['clavrx_OR_ABI-L1b-RadF-M6C01_G16_s' fileid_full '_BARBADOS-2KM-FD.level2.nc'];
                    absFN = [path2data filesep NCFN];
                    
                elseif strcmp(obj.SateSource, '2km_10min_fulldisk')
                    
                    fileID = [datestr(tnow,'yyyymmdd') HHstr MMr_str];    % this is for the fulldisk data.
                    NCFN = ['Emultic2kmNC4_goes16_' fileID '.nc'];
                    absFN = [path2data filesep NCFN];
                    
                    
                end
                
            end
            
            
            %%%%%%% Nested function 2: read subsets of data: (call different
            %%%%%%% function outside of the loop)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function datasub = read_subset_of_data(absFN, scene_center, tnow, box_width)
                if contains(absFN, 'GOES')
                    % this file contains multiple time instances.
                    datasub = read_subset_from_MPI_GOES(absFN, SateChannel, scene_center, tnow, box_width/2);
                    
                elseif contains(absFN,'clavrx')
                    % the file name contains the time stamp already
                    datasub = read_subset_from_clvarx_ABI_L1b(absFN, SateChannel, scene_center, box_width/2);
                    
                elseif contains(absFN, 'Emultic')
                    % the file name contains the time stamp already
                    datasub = read_subset_from_Emultick2kmNC4(absFN, SateChannel, scene_center, box_width/2);
                    
                end
                
            end
        end

        
        %%%% ===================   function #2     =================== %%%%
        % Purpose: construct examination region moving along a lagrangian
        % trajectory
        function lag_data_obj = construct_examination_region(obj, lagtraj, box_width)
            % the examination is constructed along the lagrangian
            % trajectory along which air flows pass the center of the SST gradient;
            % lagtraj should have the same segment numbers as those in the
            % seg data.
            
            datasub = obj.data;
            num_segs = length(datasub);
            for iseg = 1:num_segs
                % 1. interpolate the trajectory to a 10-min trajecotry
                %---- is this trajectory output as a function of time???
                % 1. find out the center of the cloud scene every 10min because
                % the resolution of the data is 10min (depends on the wind
                % speed)
                
                x0_moving = interp1(lagtraj(iseg).time, lagtraj(iseg).lon,datasub(iseg).time);
                y0_moving = interp1(lagtraj(iseg).time, lagtraj(iseg).lat,datasub(iseg).time);
                
                
                % 3. put a box around this moving center
                x_range = [x0_moving - box_width/2,  x0_moving + box_width/2];
                y_range = [y0_moving - box_width/2, y0_moving + box_width/2];
                
                lon_mask = (datasub(iseg).lon >=x_range(1)) & (datasub(iseg).lon >=x_range(2));
                lat_mask = (datasub(iseg).lat >=y_range(1)) & (datasub(iseg).lat >=y_range(2));
                
                
                % 4. take the subset of data within the box
                % do something to make sure that the size of the
                lag_data(iseg).lon{it} = datasub(iseg).lon(lon_mask);
                lag_data(iseg).lat{it} = datasub(iseg).lat(lat_mask);
                lag_data(iseg).values{it} = datasub(iseg).values(lat_mask, lon_mask);
                
            end
            
            % create dataobj for output:
            lag_data_obj = ATOMIC_GOES16();
            lag_data_obj.data = lag_data;
            obj_fields = fieldnames(obj);
            % the rest of the object properties will be the same as the
            % input object
            for i = 1:length(obj_fields)
                FN = obj_fields{i};
                if ~strcmp(FN, 'data')
                    lag_data_obj.(FN) = obj.(FN);
                end
            end
                        
            % 5. save this data;
            newFN = [obj.access(1:end-3) '_lagsubset.mat'];
            lag_data_obj.access = newFN;
            save(newFN, 'lag_data_obj');
            
        end
            
        
        %%%% ===================   function #3     =================== %%%%
        % Purpose: compute statstics for a cloud scene
        function metrics_obj = compute_cloud_metrics(obj)
            % steps: 
            %  1. identify cloud element in the scene via regionprop
            
            %  2. visualize results;
            
            %  3. compute cloud metrics
            
            
            
            %%%%%%% Nested function for sanity check on results:
            function check_regionprop_results()
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        %%%% ===================   function #4     =================== %%%%
        % Purpose: visualize cloud statistics.
        function hfig = visualize_cloud_metrics(metrics_obj)
            
        end
        
        
    end
    
    
end