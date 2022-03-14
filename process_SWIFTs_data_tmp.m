% This script will process the SWIFTs data and feed it into the wave
% coherence toolbox.

% load in SWIFTs (already in)
addpath('./func_bin');
figsvdir = './Figs/SWIFTs';
if ~exist(figsvdir)
    mkdir(figsvdir);
end

svdatadir = './processed_data';
if ~exist(svdatadir)
    mkdir(svdatadir)
end

% separate SWIFTs by legs (leg1 and leg2)
t_leg1 = [datenum(2020,1,14), datenum(2020,1,23)];
t_leg2 = [datenum(2020,1,30), datenum(2020,2,13)];


floatIDs = fieldnames(SWIFTs_ds);
NumSWIFTs = length(floatIDs);
SWIFTs_xres = 1;  % km

for ii = 1:NumSWIFTs
    FN = floatIDs{ii};
    selcrit1 = SWIFTs_ds.(FN).time>= t_leg1(1) &  SWIFTs_ds.(FN).time<=t_leg1(2);
    selcrit2 = SWIFTs_ds.(FN).time>= t_leg2(1) &  SWIFTs_ds.(FN).time<=t_leg2(2);
        
    SWF_ds_obj = ATOMIC_dataProcess(SWIFTs_ds.(FN));
    SWF_ds_obj.ATOMIC_platform = FN;
    
    SWF_leg1 = SWF_ds_obj.group_data_by_inputmask(selcrit1);
    SWF_leg2 = SWF_ds_obj.group_data_by_inputmask(selcrit2);
   
    % establish two temporary objects:
    SWF_obj_tmp1 = SWF_ds_obj;
    SWF_obj_tmp2 = SWF_obj_tmp1;
    
    % render new values for each object:
    SWF_obj_tmp1.Value = SWF_leg1;
    SWF_obj_tmp2.Value = SWF_leg2;
    
    % compute trajectory:
    SWF_obj_tmp1 = SWF_obj_tmp1.drop_bad_records;
    SWF_obj_tmp1 = SWF_obj_tmp1.compute_distance_travelled;
    
    SWF_obj_tmp2 = SWF_obj_tmp2.drop_bad_records;
    SWF_obj_tmp2 = SWF_obj_tmp2.compute_distance_travelled;
    
    mean_trajres(ii,1) = mean(diff(SWF_obj_tmp1.Value.traj),'omitnan');
    mean_trajres(ii,2) = mean(diff(SWF_obj_tmp2.Value.traj),'omitnan');
    
    % map data to an equally spaced distance axis:
    SWIFTs_remapped_obj1(ii) = SWF_obj_tmp1.map_to_distance_axis(SWIFTs_xres);   % good, I have actually did this.
    SWIFTs_remapped_obj2(ii) = SWF_obj_tmp2.map_to_distance_axis(SWIFTs_xres);   % good, I have actually did this.

    % cannot re-orient data because we don't know what the wind direction
    % is.

    
end

svdataname = ['SWIFTs_legs_ready_for_wavelet_coherence_xres' num2str(SWIFTs_xres) 'km.mat'];
save([svdatadir filesep svdataname],'SWIFTs_remapped_obj1','SWIFTs_remapped_obj2');

for ii = 1:NumSWIFTs
    % if exist wind direction parameter, then reorient data and project
    % wind;
    avail_varn = fieldnames(SWIFTs_remapped_obj1(ii).Value);
    if ismember('wind_direction', avail_varn)
        pause
        SWIFTs_remapped_obj1(ii) = SWIFTs_remapped_obj1(ii).project_wind_onto_trajectory;
        SWIFTs_remapped_obj1(ii) = SWIFTs_remapped_obj1(ii).reorient_data_sequence;
        
        SWIFTs_remapped_obj2(ii) = SWIFTs_remapped_obj1(ii).project_wind_onto_trajectory;
        SWIFTs_remapped_obj2(ii) = SWIFTs_remapped_obj1(ii).reorient_data_sequence;
    end
end


SWF_seglens={[];[]};
for j = 1:2
    eval(['inobj= SWIFTs_remapped_obj' num2str(j) ';']);
    
    for i = 1:length(inobj)
        sl = inobj(i).Value.traj(end) - inobj(i).Value.traj(1);
        SWF_seglens{j} = [SWF_seglens{j}, sl];
    end
end


figure(7);
hold on;
for j = 1:2
    bar(SWF_seglens{j},'FaceAlpha',0.5)
end

% call into the wavelet coherence toolbox
% ideally, the wave direction should be rather similar with the wind
% direction if the measured waves are wind waves. 
% how to understand the wind speed measured by drifters and other
% instruments???
varnlist = { 'wind_speed', 'wind_direction','u_algtraj'};
svfig = false;

NumSWIFTs = length(SWIFTs_remapped_obj1);
for iv = 1:length(varnlist)
    wind_varn = varnlist{iv};
    
    clear SWF_wtc_stat SWF_outCOI_stat SWF_inCOI_stat
    for legnum=1:2
        jj=0;
        %legnum = 1;
        for i = 1:NumSWIFTs
            
            eval(['inobj = SWIFTs_remapped_obj' num2str(legnum) '(i);']);     % how would the sampling resolution and remapping resolution change the results?
            inobj.Value
            for iseg =1:length(inobj.Value)
                jj = jj+1;
                
                ts1 =inobj.Value(iseg).sea_water_temperature;
                try
                    ts2 = inobj.Value(iseg).(wind_varn);
                    traj = inobj.Value(iseg).distance_axis;
                    time = inobj.Value(iseg).time;
                    
                    validID = ~isnan(time) & ~isnan(ts2) & ~isnan(ts1);
                    ts1 = ts1(validID);
                    ts2= ts2(validID);
                    traj = traj(validID);
                    time = time(validID);
                    
                    
                    ATOMIC_platform = inobj.ATOMIC_platform;
                    
                    labelstr.ts1 = 'SST (^{\circ}C)';
                    labelstr.ts2 = {'wind speed (m/s)'}; %'along-traj.';
                    
                    if (max(traj) - min(traj))>80
                        
                        %%% 0. establish the wavelet object for each glider:
                        SWF_wtcobj = ATOMIC_Wavelet(ts1, ts2, traj, time, ATOMIC_platform, SWIFTs_xres);
                        hfig0 = SWF_wtcobj.plot_data(labelstr);
                        if svfig
                            xc_savefig(hfig0,figsvdir, [ATOMIC_platform '_leg' num2str(legnum) ...
                                '_SST_' wind_varn 'record.jpg'], [0, 0, 10 8]);
                        end
                        
                        
                        %%% 1. apply box-cox transform on the data to have normal distribution;
                        [SWF_wtcobj_trans, hfig1] = SWF_wtcobj.transform_to_normal_distribution(labelstr);
                        if svfig
                            xc_savefig(hfig1,figsvdir, [ATOMIC_platform '_leg' num2str(legnum) ...
                                '_SST_' wind_varn 'record_boxcox_transformed.jpg'], [0, 0, 10 8]);
                        end
                        %
                        
                        %%% 2. apply the wave coherence toolbox:
                        % [hf'part-' num2str(cnt)igs, wtc_stat(iseg)] = WG247_wtcobj.wavelet_coherence_toolbox;
                        
                        [hfigs, SWF_wtc_stat(jj, legnum), ...
                            SWF_outCOI_stat(jj,legnum), SWF_inCOI_stat(jj, legnum)] = SWF_wtcobj_trans.wavelet_coherence_toolbox;
                        
                        figure(hfigs.Number);
                        title([ATOMIC_platform ': ' datestr(time(1)) '~' datestr(time(end))]);
                        
                        if svfig
                            xc_savefig(gcf,figsvdir, [ATOMIC_platform '_leg' num2str(legnum) ...
                                '_coherence__boxcox_transformed_SST_and_' wind_varn '.jpg'], [0, 0, 10 8]);
                        end
                        pause (0.5);
                        % close all;
                        close all;
                    end
                    
                end
            end
            
        end
    end
    
    % save output for plotting purpose;
    try
        svdataname = ['SWIFTs_wavecoherence_stat_SST_and_' wind_varn '_xres' num2str(SWIFTs_xres) 'km.mat'];
        save([svdatadir filesep svdataname],'SWF_wtc_stat','SWF_outCOI_stat','SWF_inCOI_stat');
    end
    
end
%

