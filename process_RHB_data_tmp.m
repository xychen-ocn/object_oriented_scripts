% This script is used to process Ron Brown data; 

% load in the RHB datasets:
addpath('./func_bin');
figsvdir = './Figs/RHB';
if ~exist(figsvdir)
    mkdir(figsvdir);
end

svdatadir = './processed_data';
if ~exist(svdatadir)
    mkdir(svdatadir)
end


% extract data on the selected days:
DOIs ={'Jan08','Jan09','Jan10',...
       'Jan18',...
       'Jan22','Jan23','Jan24', ...       
       'Feb03','Feb04', ...
       'Feb06','Feb07','Feb08', ...
       'Feb10','Feb11','Feb12'};
   
seg_times={{'Jan08','Jan11'}, {'Jan18', 'Jan19'}, {'Jan22','Jan25'}, ...
          {'Feb03','Feb05'}, {'Feb06','Feb09'},{'Feb10','Feb13'}};

% what is RHB_ds?? where did it come from?
      
nseg = length(seg_times);
for iseg = 1:nseg
    seg_start_t = datenum(['2020' seg_times{iseg}{1}],'yyyymmmdd');
    seg_end_t   = datenum(['2020' seg_times{iseg}{2}],'yyyymmmdd');
    seg_mask = (RHB_ds.time-4/24>=seg_start_t & RHB_ds.time-4/24<=seg_end_t);       % RHB_ds.time is UTCtime. 
    
    RHB_ds_obj = ATOMIC_dataProcess(RHB_ds);
    RHB_ds_obj.ATOMIC_platform = 'RHB';
    RHB_grped_segs(iseg) = RHB_ds_obj.group_data_by_inputmask(seg_mask);
    % 0. drop bad records:
    RHB_obj_tmp = RHB_ds_obj;
    RHB_obj_tmp.Value = RHB_grped_segs(iseg);
    RHB_obj_tmp_good(iseg) = RHB_obj_tmp.drop_bad_records;   
    
    % 1. compute distance travelled :
    RHB_obj_tmp_good(iseg) = RHB_obj_tmp_good(iseg).compute_distance_travelled;

end
for iseg = 1:nseg
    mean_RHBxres(iseg)= mean(diff(RHB_obj_tmp_good(iseg).Value.traj),'omitnan');
end

% 2. select straight traj elements
crits.ddir = 15;
crits.misang=60;     % maximum deviation angle from the along-wind axis. 
crits.length= 100;   % units: km
sstname='tsea';
RHB_xres=2;          % units: km

for iseg = 1:nseg
    % remap data to the distance axis first:
    %RHB_obj_remapped = RHB_obj_tmp_good(iseg);
    RHB_obj_remapped = RHB_obj_tmp_good(iseg).map_to_distance_axis(RHB_xres);
    
    % select straight segments:
    RHB_straight_segs(iseg) = RHB_obj_remapped.select_straight_traj_segments(crits, sstname);
    % save figure:
%     figname =['RHB_seg' num2str(iseg,'%2.2i') '_' seg_times{iseg}{1} '_' seg_times{iseg}{2} ...
%               '_breaking_down_data_crit-ddir' num2str(crits.ddir) ...
%               '_crit-misang' num2str(crits.misang) '_crit-length' num2str(crits.length) '.jpg'];
% 
%     xc_savefig(gcf,figsvdir,figname,[0 0 12 10]);
end


cnt = 0;
for iseg = 1:nseg
    
    if ~isempty(RHB_straight_segs(iseg).along_wind)
        cnt=cnt+1;
        % 3. project wind on to the trajectory
        RHB_alongwind_segobj(cnt) = ATOMIC_dataProcess(RHB_straight_segs(iseg).along_wind);
        RHB_alongwind_segobj(cnt) = RHB_alongwind_segobj(cnt).project_wind_onto_trajectory;
        
        % 3.b compute Ro number based on the along trajectory vorticity.
        RHB_alongwind_segobj(cnt) = RHB_alongwind_segobj(cnt).compute_along_traj_current_relative_vorticity_and_Ro;
        % 4. reoient data sequence
        RHB_alongwind_segobj(cnt)= RHB_alongwind_segobj(cnt).reorient_data_sequence;
    end
   % close all
    %pause
end

% check the average segment length;  300 ~550km, 2segs less than 150m  
RHB_seglens = [];
for i = 1:length(RHB_alongwind_segobj)
    RHB_seglens = [RHB_seglens, [RHB_alongwind_segobj(i).Value.segment_length]];
end

figure
bar(RHB_seglens)



% break down data into records with the same heading, and re-orient data;
% save the data at this step;
% save the RHB processed data now:
svdataname = ['RHB_along_wind_segments_ready_for_wavelet_coherence_updated_Ro_' date '.mat'];   %_crit80m
save([svdatadir filesep svdataname],'RHB_alongwind_segobj','RHB_straight_segs','RHB_seglens');



load([svdatadir filesep svdataname]);

%% put data into the wavelet coherence toolbox.
varnlist = { 'wspd_10N', 'wdir','u_algtraj'};
svfig = false;
for iv = 1:length(varnlist)
    wind_varn = varnlist{iv};
    
    clear RHB_wtc_stat RHB_inCOI_stat RHB_outCOI_stat
    jj=0;
    for i = 1:length(RHB_alongwind_segobj)
        inobj = RHB_alongwind_segobj(i);     % how would the sampling resolution and remapping resolution change the results?
        
        for iseg =1:length(inobj.Value)
            jj = jj+1;
            %validID = ~isnan(ts1);
            ts1 =inobj.Value(iseg).tsea;
            %ts1= ts1
            %ts2 = WG247_alongwind_segobj_remapped.Value(iseg).wind_speed;
            ts2 = inobj.Value(iseg).(wind_varn);
            traj = inobj.Value(iseg).distance_axis;
            time = inobj.Value(iseg).time;
            
            ts_pair(jj).ts1 = ts1;
            ts_pair(jj).ts2 = ts2;
            ts_pair(jj).traj = traj;
            ts_pair(jj).time = time;
        
            ATOMIC_platform = 'RHB';
            
            labelstr.ts1 = 'SST (^{\circ}C)';
            labelstr.ts2 = {'along-traj.','wind speed (m/s)'}; %'along-traj.';
            
            %%% 0. establish the wavelet object for each glider:
            RHB_wtcobj = ATOMIC_Wavelet(ts1, ts2, traj, time, ATOMIC_platform, RHB_xres);
            
            if svfig
                hfig0 = RHB_wtcobj.plot_data(labelstr);
                xc_savefig(hfig0,figsvdir, [ATOMIC_platform '_seg' num2str(i, '%2.2i') ...
                    'part-' num2str(iseg) '_SST_' wind_varn 'record.jpg'], [0, 0, 10 8]);
            end
            
            
            %%% 1. apply box-cox transform on the data to have normal distribution;
            [RHB_wtcobj_trans, hfig1] = RHB_wtcobj.transform_to_normal_distribution(labelstr);
            
            if svfig
                xc_savefig(hfig1,figsvdir, [ATOMIC_platform '_seg' num2str(i, '%2.2i') ...
                    'part-' num2str(iseg) '_SST_' wind_varn 'record_boxcox_transformed.jpg'], [0, 0, 10 8]);
            end
                    %
            
            %%% 2. apply the wave coherence toolbox:
            % [hf'part-' num2str(cnt)igs, wtc_stat(iseg)] = WG247_wtcobj.wavelet_coherence_toolbox;
            
            [hfigs, RHB_wtc_stat(jj), RHB_outCOI_stat(jj), RHB_inCOI_stat(jj)] = RHB_wtcobj_trans.wavelet_coherence_toolbox;
            
            figure(hfigs.Number);
            title([ATOMIC_platform ': ' datestr(time(1)) '~' datestr(time(end))]);    %UTC time;
            
            if svfig
                xc_savefig(gcf,figsvdir, [ATOMIC_platform '_seg' num2str(i, '%2.2i') ...
                    'part-' num2str(iseg) '_coherence__boxcox_transformed_SST_and_' wind_varn '.jpg'], [0, 0, 10 8]);
            end
            pause (0.5);
            close all;
            
            
            
            %         %% make use of information obtained from wtc_stat:
            %         [xx, yy] = meshgrid(traj, log2(RHB_wtc_stat(1).period));
            %         Rsq = RHB_wtc_stat(1).coherence_squared;
            %         Rsq_locmax = imregionalmax(Rsq);
            %         sig95_locs = RHB_wtc_stat(1).sig95 >=1;  % the definition of sig95 was wtcsig = Rsq./wtcsig in the code..
            %         % remove locs within the COI:
            %         outcoi_locs = Rsq > RHB_wtc_stat(1).coi;  % wavelength as a function of physical location.
            %
            %         % use inpolygon:  create a polygon area for the COI;
            %         % if the sig95_locs are within polygon, then collect and label the
            %         % data separately. (do not discard it completely..)
            %
            %
            %         % obtain averaged phase and scale within the regions of interests:
            %         % xv, yv represent the contour (polygon) of the 95% significant
            %         % region. xx, yy is the two axis I had.
            %
            %         % use regionprop to find the area with 95% significance:
            %         CCr = bwconncomp(sig95_locs,8);
            %         Lr = labelmatrix(CCr);
            %         stats=regionprops(Lr, sig95_locs,'PixelIdxList','Image');
            %         num_sig95_patches = length(stats);
            %
            %         % use pixel index to compute the averaged phase and ...
            %         xcoi = [traj([1 1])-xres*.5, traj, traj([end end])+xres*.5];
            %         ycoi = log2([RHB_wtc_stat(1).period([end 1]) ...
            %             RHB_wtc_stat(1).coi ...
            %             RHB_wtc_stat(1).period([1 end])]);
            %
            %         wvlen_scale = RHB_wtc_stat(1).period.*xres;
            %         [~, WVLEN] = meshgrid(traj, wvlen_scale);
            %
            %         cnt1=0;
            %         cnt2=0;
            %         for ii = 1:num_sig95_patches
            %             idx = stats(ii).PixelIdxList;
            %             xsel = xx(idx);
            %             ysel = yy(idx);
            %
            %             [inCOI, onCOI] = inpolygon(xsel, ysel, xcoi, ycoi);
            %
            %             cond = ~inCOI & ~onCOI ;
            %             num_pixel_outCOI = length(idx(cond));
            %
            %             if num_pixel_outCOI>0.8*length(cond)
            %                 % significant 95% location is outside of the COI;
            %                 cnt1 =cnt1+1;
            %                 notCOI.ave_phase(cnt1) = mean(RHB_wtc_stat(1).phase_angle(idx));
            %                 notCOI.ave_wvlen(cnt1) = mean(WVLEN(idx));
            %             else
            %                 cnt2=cnt2+1;
            %
            %                 COI_area(cnt2).ave_phase = mean(RHB_wtc_stat(1).phase_angle(idx));
            %                 COI_area(cnt2).ave_wvlen = mean(WVLEN(idx));
            %                 COI_area(cnt2).fraction = num_pixel_outCOI./length(cond);
            %                 figure(1);
            %                 plot(xsel, ysel,'xk');
            %             end
            %         end
            %
            %         figure;
            %         plot([notCOI.ave_wvlen], [notCOI.ave_phase],'.b');
            %
        end
        
    end
    
    timestamp = datestr(now, 'yyyy-mmm-dd');
    % save output for plotting purpose;
    svdataname = ['RHB_wavecoherence_stat_SST_and_' wind_varn '_xres' num2str(RHB_xres) 'km_' timestamp '.mat'];
    save([svdatadir filesep svdataname],'RHB_wtc_stat','RHB_outCOI_stat','RHB_inCOI_stat','ts_pair');
end

%
figure(21);clf;
scatter([outCOI_stat.ave_wvlen], [outCOI_stat.ave_phase].*180/pi, 30, [outCOI_stat.ave_cohsq], 'filled');
hold on;
scatter([inCOI_stat.ave_wvlen], [inCOI_stat.ave_phase].*180/pi,'marker','x');



% I may need to write a README.md to describe how to use these
% object-oriented scripts. 
% visualize the surface current:
% visobj = ATOMIC_dataVisual();
% for iseg = 1:5
%     visobj.data = RHB_alongwind_segobj(iseg).Value;
%     visobj.fignum = iseg;
%     visobj.plot_current_along_track;
% end
