% This script will be used to handle the data processing for the
% wavegliders.
% 

% assume that the wave glider data has been loaded already;
% WaveGlider_ds;

figsvdir = './Figs/WGs';
if ~exist(figsvdir)
    mkdir(figsvdir);
end

svdatadir = './processed_data';
if ~exist(svdatadir)
    mkdir(svdatadir)
end


% 1. compute the heading of wave gliders: (compute from the distance
% traveled)
WGIDs = fieldnames(WaveGlider_ds);
NumWGs = length(WGIDs);

for i = 1:NumWGs
    FN = WGIDs{i};
    obj_tmp = ATOMIC_dataProcess(WaveGlider_ds.(FN));    % temporary object;
    if i==1
        obj_tmp = obj_tmp.drop_bad_records;
        for j = 1:length(obj_tmp.Value)
            obj_part = ATOMIC_dataProcess(obj_tmp.Value(j));
            WGs_obj.(FN)(j) = obj_part.compute_distance_travelled;
        end
    else
        
        WGs_obj.(FN) = obj_tmp.compute_distance_travelled;
    end
end

% check the headings of each wave gliders:
figure(gcf)
subplot(211)
plot(WGs_obj.WG245(1).Value.traj, WGs_obj.WG245(1).Value.trajdir_cart);
hold on;
plot(WGs_obj.WG245(2).Value.traj, WGs_obj.WG245(2).Value.trajdir_cart);

subplot(212)
plot(WGs_obj.WG247.Value.traj, WGs_obj.WG247.Value.trajdir_cart,'.-k');


% figure
% plot(WGs_obj.WG247.Value.traj, WGs_obj.WG247.Value.wind_direction);


% break the WG records apart according to the following criteria:
% 1. direction change of the trajectory less than 10° (need a straight line)
% 2. the direction between the trajectory and the wind is not within
% 3. length of the records > 100km or 150km (after segment is selected.)
% (90+-10°) (a.k.a wind is mainly along the axis of the trajectory.)


% 1.
crits.ddir = 10;
crits.misang=60;     % maximum deviation angle from the along-wind axis. 
crits.length= 100;   % units: km


WGs_obj2 = ATOMIC_dataProcess(WGs_obj.WG247.Value);
WGs_obj2.ATOMIC_platform='WG247';
WG247_segs = WGs_obj2.select_straight_traj_segments(crits,'sea_water_temperature');
figname =[WGs_obj2.ATOMIC_platform '_breaking_down_data_crit-ddir' num2str(crits.ddir) ...
     '_crit-misang' num2str(crits.misang) '_crit-length' num2str(crits.length) '.jpg'];
xc_savefig(gcf,figsvdir,figname, [0 0 12 10]);

% need to determine whether or not the segment needs to be re-oriented
% (this will be done using the dataProcess class;)

WG247_alongwind_segobj = WGs_obj2;
WG247_alongwind_segobj.Value = WG247_segs.along_wind;
WG247_alongwind_segobj = WG247_alongwind_segobj.project_wind_onto_trajectory;
WG247_alongwind_segobj = WG247_alongwind_segobj.reorient_data_sequence;
xres = 1;   % units: km
WG247_alongwind_segobj_remapped = WG247_alongwind_segobj.map_to_distance_axis(xres);

svdataname = [WG247_alongwind_segobj.ATOMIC_platform '_along_wind_segments_ready_for_wavelet_coherence.mat'];
save([svdatadir filesep svdataname],'WG247_alongwind_segobj','WG247_segs','WG247_alongwind_segobj_remapped');


%% WG245:
for ik = 1:2
    WGs_obj2 = ATOMIC_dataProcess(WGs_obj.WG245(ik).Value);
    WGs_obj2.ATOMIC_platform=['WG245-part' num2str(ik)];
    
    %     WG245_segs(ik) = WGs_obj2.select_straight_traj_segments(crits,'sea_water_temperature');
    %     figname =[WGs_obj2.ATOMIC_platform '_breaking_down_data_crit-ddir' num2str(crits.ddir) ...
    %      '_crit-misang' num2str(crits.misang) '_crit-length' num2str(crits.length) '.jpg'];
    %     xc_savefig(gcf,figsvdir, figname, [0 0 12 10]);
    %
    %  need to determine whether or not the segment needs to be re-oriented
    % (this will be done using the dataProcess class;)
    if ~isempty(WG245_segs(ik).along_wind)
        WG245_alongwind_segobj = WGs_obj2;
        WG245_alongwind_segobj.Value = WG245_segs(ik).along_wind;
        WG245_alongwind_segobj = WG245_alongwind_segobj.project_wind_onto_trajectory;
        WG245_alongwind_segobj = WG245_alongwind_segobj.reorient_data_sequence;
        % interpolate data to a equally spaced distance axis;
        xres = 1;   % units: km
        WG245_alongwind_segobj_remapped = WG245_alongwind_segobj.map_to_distance_axis(xres);
        
        % save the processed data from above:
        svdataname = [WG245_alongwind_segobj.ATOMIC_platform '_along_wind_segments_ready_for_wavelet_coherence.mat'];
        save([svdatadir filesep svdataname],'WG245_alongwind_segobj','WG245_segs','WG245_alongwind_segobj_remapped');
        
        
    end
    % pause
end


%% put these segment into the wavelet coherence scripts:


inobj = WG245_alongwind_segobj_remapped;     % how would the sampling resolution and remapping resolution change the results?

wind_varn = 'wind_direction';
for iseg = 1:length(inobj.Value)
    %validID = ~isnan(ts1);
    ts1 =inobj.Value(iseg).sea_water_temperature;
    %ts1= ts1
    %ts2 = WG247_alongwind_segobj_remapped.Value(iseg).wind_speed;
    ts2 = inobj.Value(iseg).(wind_varn);
    traj = inobj.Value(iseg).distance_axis;
    time = inobj.Value(iseg).time;
    

    ATOMIC_platform = inobj.ATOMIC_platform;
    
    labelstr.ts1 = 'SST (^{\circ}C)';
    labelstr.ts2 = {'along-traj.','wind speed (m/s)'}; %'along-traj.';

    %%% 0. establish the wavelet object for each glider:
    WG_wtcobj = ATOMIC_Wavelet(ts1, ts2, traj, time, ATOMIC_platform, xres);
    hfig0 = WG_wtcobj.plot_data(labelstr);
    xc_savefig(hfig0,figsvdir, [ATOMIC_platform '_seg' num2str(iseg, '%2.2i') ...
        '_SST_' wind_varn '_record.jpg'], [0, 0, 10 8]);

    
    %%% 1. apply box-cox transform on the data to have normal distribution;
    [WG_wtcobj_trans, hfig1] = WG_wtcobj.transform_to_normal_distribution(labelstr);
    
    xc_savefig(hfig1,figsvdir, [ATOMIC_platform '_seg' num2str(iseg, '%2.2i') ...
        '_SST'  wind_varn  '_record_transformed_boxcox.jpg'], [0, 0, 10 8]);
    
    
    %%% 2. apply the wave coherence toolbox:
       % [hfigs, wtc_stat(iseg)] = WG247_wtcobj.wavelet_coherence_toolbox;
    tmpN = strsplit(ATOMIC_platform,'-');
    PN = [];
    for k = 1:length(tmpN)
        PN = [PN tmpN{k}];
    end
    
    [hfigs, WG_wtc_stat.(PN)(iseg), WG_outCOI_stat.(PN)(iseg), WG_inCOI_stat.(PN)(iseg)] = WG_wtcobj_trans.wavelet_coherence_toolbox;

    figure(hfigs.Number);
    title([ATOMIC_platform ': ' datestr(time(1)) '~' datestr(time(end))]);
    
    
    xc_savefig(gcf,figsvdir, [ATOMIC_platform '_seg' num2str(iseg, '%2.2i') ...
        '_coherence_transformed_boxcox_SST_and_' wind_varn '.jpg'], [0, 0, 10 8]);
    
    pause(0.5)
    close all;
    
    
    
    
    %% make use of information obtained from wtc_stat:
    
    

    
end


% save output for plotting purpose;
svdataname = ['WGs_wavecoherence_stat_SST_and_' wind_varn '.mat'];
save([svdatadir filesep svdataname],'WG_wtc_stat','WG_outCOI_stat','WG_inCOI_stat');

%% note: several thing to be tested:
% 1. what resolution should be used for the interpolation? will it affect
% the coherence results?

% 2. total wind speed versus wind speed projected on to the SST gradient.
%    verus wind speed normal to the trajectory;

% 3. ts2 = wind direction.

% 4. are there common clusters in the scatter plots??
%    - find local maximum in the coherence
%    - conditional selection: select local coherence maximum that is higher
%    than the sig95 threshold.
%    - get the averaged phase angle and wavelength from the region with
%    high coherence
%    - plot it on a scatter plots. 
