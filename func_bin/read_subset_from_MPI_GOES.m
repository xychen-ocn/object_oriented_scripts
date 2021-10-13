function datasub = read_subset_from_MPI_GOES(NCFN, channel, scene_center, tnow, radius)
%% purpose: take subset from GOES16 satellite data postprocessed by Hauke(MPI).
%% input: filename, channel, scene_center, raidus
%% output: subsetted data (struct)

% purpose: carve out a spatial box center at (x0, y0) in a square with a side of
% 2xradius; take the data in center at t0, with a give time window.
% (t-window/2: t+window/2)
% input: data  (a matlab structure containing coordinate, time and var)
% output: datasub (subsampled data in space and time)
%

x0 = scene_center(1);
y0 = scene_center(2);

lon_range = [x0-radius, x0+radius];
lat_range = [y0-radius, y0+radius];
%time_range = tnow;


if strcmpi(channel, 'VIS')
    varn='C02';
else
    varn='C13';
end

% read in data:
try
    [data, errmsg] = read_ncfile_to_struct(NCFN.default);
catch  ME
    disp(ME.identifier)
    disp(ME.message)
    if cnotains(ME.message, '-90')
        % file not found, use another data source;
        disp('using another source file from 2km_10min folder');
        default_flag = false;
        
    end
    %             if strcmp(channel, 'IR')
    %                 % try the other source:
    %                 NCFN0 =  [datadir filesep datestr(datehere,format1) filesep ...
    %                     'GOES_02_8N-18N-62W-50W_' datestr(datehere,format2) '.nc'];
    %                 [data, errmsg] = read_ncfile_to_struct(NCFN0);
    %             end
end

for iv = 1:length(errmsg)
    if ~isempty(errmsg{iv})
        disp(['warning: ' errmsg{iv}]);
    end
end


[~, lon_stid] = min(abs(data.lon - lon_range(1)));
[~, lon_edid] = min(abs(data.lon - lon_range(2)));
[~, lat_stid] = min(abs(data.lat - lat_range(1)));
[~, lat_edid] = min(abs(data.lat - lat_range(2)));

[~, time_stid] = min(abs(data.time - tnow));


% another way is to use mask;
% lon_mask = (data.lon>=lon_range(1)) .* (data.lon <=lon_range(2));
% lat_mask = (data.lat>=lat_range(1)) .* (data.lat <=lat_range(2));
% time_mask = (data.time>=time_range(1)) .* (data.time <=time_range(2));

% use the index to extract the subset;
datasub.lon = data.lon(lon_stid:lon_edid);
datasub.lat = data.lat(lat_stid:-1:lat_edid);
datasub.time = data.time(time_stid);

% take a subset of the C02/C13 data:
% read from the NCFN again:
% nt = time_edid - time_stid +1;
% disp(['--> extracting data from '  datestr(data.time(time_stid)) ' to ' ...
%       datestr(data.time(time_edid))]);
% disp(['--> num of cloud scenes: ' num2str(ceil(nt/tinc)) ]); 

tmp = ncread(NCFN, varn, [lon_stid, lat_edid, time_stid],[lon_edid-lon_stid+1, abs(lat_edid-lat_stid)+1, 1]);
datasub.values = permute(tmp(:,end:-1:1),[2,1]);   % transpose/permute to NLAT X NLONend







end