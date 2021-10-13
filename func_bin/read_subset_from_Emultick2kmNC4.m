function read_subset_from_Emultick2kmNC4(NCFN, channel, scene_center, radius)
%% purpose: take subset from clvarx_ABI_L1b data on the aeris site.
%% input: filename, channel, scene_center, raidus
%% output: subsetted data (struct)

x0 = scene_center(1);
y0 = scene_center(2);

% read in the netCDF file:
X2km = ncread(NCFN, 'X2km');
Y2km = ncread(NCFN, 'Y2km');

lat_proj_origin = ncreadatt(NCFN,'geos','latitude_of_projection_origin');
lon_proj_origin = ncreadatt(NCFN,'geos','longitude_of_projection_origin');

lat_pc = (Y2km-0)./111E3 + lat_proj_origin;
lon_pc = (X2km-0)./(111E3.*cosd(lat_pc)) + lon_proj_origin;


if strcmp(channel, 'VIS')
    %varn = 'refl_0_65um_nom';
    varn = 'VIS_006';
    
else
    %varn = 'temp_10_4um_nom';
    varn = 'IR_103';
    
end
dataval = ncread(NCFN, varn);    % nx, ny

%% 3. extract a subset of data points:
lon_subrange = [x0-radius, x0+radius];
lat_subrange = [y0-radius, y0+radius];

lon_mask = (lon_pc>=lon_subrange(1)) & (lon_pc <=lon_subrange(2));
lat_mask = (lat_pc>=lat_subrange(1)) & (lat_pc <=lat_subrange(2));


% use the index to extract the subset;
datasub.lon = lon_pc(lon_mask);
datasub.lat = lat_pc(lat_mask);
datasub.values = dataval(lon_mask, lat_mask)';

end