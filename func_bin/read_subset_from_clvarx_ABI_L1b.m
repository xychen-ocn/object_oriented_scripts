function datasub = read_subset_from_clvarx_ABI_L1b(NCFN, channel, scene_center, radius)
%% purpose: take subset from clvarx_ABI_L1b data on the aeris site.
%% input: filename, channel, scene_center, raidus
%% output: subsetted data (struct)


x0 = scene_center(1);
y0 = scene_center(2);

% read in the netCDF file:
lon_pc = ncread(NCFN, 'longitude_pc');
lat_pc = ncread(NCFN, 'latitude_pc');

% construct a new equally-spaced grid to regrid the L1b data onto:
lon_range = [-62, -50];
lat_range = [8, 18];
lonres = 2E3/(111E3*cosd(10));
latres = 2E3/111E3;
lon_grd = lon_range(1):lonres:lon_range(2);
lat_grd = lat_range(1):latres:lat_range(2);
[LON_GRD, LAT_GRD]= meshgrid(lon_grd, lat_grd);

if strcmp(channel, 'VIS')
    varn = 'refl_0_65um_nom';
    
else
    varn = 'temp_10_4um_nom';
    
end
dataval = ncread(NCFN, varn);

valid = ~isnan(dataval);  % change to dataval.
F= scatteredInterpolant(lon_pc(valid), lat_pc(valid), dataval(valid),...
    'linear','none');
dataval_gridded = F(LON_GRD, LAT_GRD);

%% 3. extract a subset of data points:
lon_subrange = [x0-radius, x0+radius];
lat_subrange = [y0-radius, y0+radius];

lon_mask = (lon_grd>=lon_subrange(1)) & (lon_grd <=lon_subrange(2));
lat_mask = (lat_grd>=lat_subrange(1)) & (lat_grd <=lat_subrange(2));


% use the index to extract the subset;
datasub.lon = lon_grd(lon_mask);
datasub.lat = lat_grd(lat_mask);
datasub.values = dataval_gridded(lat_mask, lon_mask);

end