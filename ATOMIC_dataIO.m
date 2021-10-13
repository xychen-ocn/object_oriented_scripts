% ======================================================================= %
% class: ATOMIC_dataIO
% purpose: handle reading and loading ATOMIC netCDF data (local & remote)
% Author: Xuanyu Chen 
% Date: v0 @ 09/28/2021
% ======================================================================= %

classdef ATOMIC_dataIO
    properties
        DataPath 
        DataName 
        DataLoc 
        saveflag (1,1) logical 
        DataTag 
        svdatadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
    end
    
    methods 
        
        %%% ======== initialize class or obj ========== %%%
        function obj = ATOMIC_dataIO(val)
            
            if nargin> 0
                PropNames = fieldnames(val);
                for i = 1:length(PropNames)
                    PN = PropNames{i};
                    obj.(PN) = val.(PN);
                end
                
            else
                obj.DataPath = 'none';
                obj.DataName = 'none';
                obj.DataLoc = 'none';
                obj.saveflag = 'false';
                obj.DataTag = 'none';
            end
        end
        %%% ======= ncdisp ====== %%%
        function view_meta_data(obj)
            if ~contains([obj.DataName],'*')
                
                ncdisp([obj.DataPath filesep obj.DataName]);
            else
                
                filenames = get_series_of_filenames(obj);
                FN = filenames{1};    % just open one file to look at the meta data. 
                ncdisp([obj.DataPath filesep FN]);
            end
        end
        
        
        %%% ======= this is a subfunction for read_netCDF_to_metdata ==== %
        % this function get the data filename that is represented by the
        % wildcard ...
        function filenames =get_series_of_filenames(obj)
            
            % read the filenames from the data directory
            if strcmp(obj.DataLoc,'local')
                files = dir([obj.DataPath filesep obj.DataName]);
                filenames = {files.name};
                
            elseif strcmp(obj.DataLoc, 'opendap')
                % reconstruct opendap datapath:
                url = [obj.DataPath filesep 'catalog.html'];
                url  = replace(url, 'dodsC','catalog');
                
                webdata = webread(url);
                
                % get part of the data filename (before the wild card)
                FNinfo = strsplit(obj.DataName, '*');
                expression = [FNinfo{1} '\d+' FNinfo{2}];               % any digits of a numeric number
                filenames = unique(regexp(webdata, expression,'match'));      % return unique filenames
                
            end
            
        end
        
        %%% ======== read ======== %%%
        function data_structure = read_netCDF_to_matdata(obj)
            % This will be a function that read netCDF into a local matlab
            % file.
            
            %%%  check if there is wild card character in the object
            %%% (string)
            if contains([obj.DataName],'*')
                disp('--> reading a bunch of files into a structure..');
                % read a bunch of files:
                filenames = get_series_of_filenames(obj);
                
                for i = 1:length(filenames)
                    FN = filenames{i};
                    absFN = [obj.DataPath filesep FN];
                    %ncdisp(absFN);
                    
                    
                    tmp = strsplit(FN, '_');
                    % manually construct for now:
                    if strcmp(obj.DataTag, 'P3_IR')
                        dateIDstr = tmp{end-1};
                        if length(dateIDstr) ==8
                            format = 'yyyymmdd';
                        end
                        dateIDnum = datenum(dateIDstr, format);
                        DN = datestr(dateIDnum, 'mmmdd');
                        
                    elseif strcmp(obj.DataTag, 'SWIFTs') | strcmp(obj.DataTag, 'WaveGlider')
                        DN = tmp{3};
                        
                    end
                    
                    data_structure.(DN) = read_ncfile_to_struct(absFN);
                    

                end
                
            else
                disp('--> reading just ONE single file into a structure..');
                % read 1 single files:
                absFN = [obj.DataPath filesep obj.DataName];
                data_structure = read_ncfile_to_struct(absFN);
            end
            
            % save data if requested.
            if obj.saveflag
                disp('--> saving data :D');
                tmp = strsplit(obj.DataName, '_');
                versionstr = tmp{end}(1:end-3);
                if strcmp(obj.DataLoc, 'local') 
                    svdatadir = obj.DataPath;
                
                else
                    svdatadir = obj.svdatadir;
                end
                
                abs_matFN = [svdatadir filesep 'EUREC4A_ATOMIC_' obj.DataTag '_AllinOne_' versionstr '.mat'];

                eval([obj.DataTag '= data_structure;']);
                eval(['save(abs_matFN, ''' obj.DataTag ''', ''-v7.3'')']);
            end
            
        end
        
        
        %%% ========  write ======== %%%
        function data = load_local_matdata(obj)
           % This function read the locally stored matlab data 
           
           %%%%  check if data exist:
           absFN = [obj.DataPath filesep obj.DataName];
           disp(absFN)
           if exist(absFN,'file')~=0
               data = load(absFN);
               
           else
               disp('Warning: no such file! Please check input file name!');
               return
           end
            
            
        end
        
        

        

        
    end
    
    
end