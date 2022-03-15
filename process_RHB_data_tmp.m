% This script is used to process Ron Brown data; 
clear all; clc; close all;

% what is RHB_ds?? where did it come from? (I think I loaded the dataset in
% the main program.)
load('RHB_dataset.mat');

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
crits.ddir = 15;% 15
crits.misang=45;     % 60; maximum deviation angle from the along-wind axis. 
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


%% get RHB alongwind segments: (wind cross SST gradient)
cnt = 0;
for iseg = 1:nseg
    
    if ~isempty(RHB_straight_segs(iseg).along_wind)
        cnt=cnt+1;
        % 3. project wind on to the trajectory
        % initialize object:
        RHB_alongwind_segobj(cnt) = ATOMIC_dataProcess(RHB_straight_segs(iseg).along_wind);
        RHB_alongwind_segobj(cnt).ATOMIC_platform = 'RHB';
        
        % project two types of wind:
        % I seem to have problems with the sign of the wind in the latest
        % version that I created. 
        % What I need: after projection, the wind component in the downwind
        % direction should be postive.

        RHB_alongwind_segobj(cnt) = RHB_alongwind_segobj(cnt).project_wind_onto_trajectory('wspd_varn','wspd');     % I seem to have problems with the sign of the wind in the latest version that I created.
        RHB_alongwind_segobj(cnt) = RHB_alongwind_segobj(cnt).project_wind_onto_trajectory('wspd_varn','wspd_10N');
        
        % 3.b compute Ro number based on the along trajectory vorticity.
        RHB_alongwind_segobj(cnt) = RHB_alongwind_segobj(cnt).compute_along_traj_current_relative_vorticity_and_Ro;
        
        % 4. reoient data sequence (should this be the last operation to the data? -> Yes,because most of the varialbes will be reoriented if needed.)
        RHB_alongwind_segobj(cnt)= RHB_alongwind_segobj(cnt).reorient_data_sequence;

        % after re-orientation, compute the gradient in the downwind
        % direction.
        % compute SST gradient and laplacian using (central difference):
        RHB_alongwind_segobj(cnt)= RHB_alongwind_segobj(cnt).compute_grad_Laplacian('tsea');
        RHB_alongwind_segobj(cnt)= RHB_alongwind_segobj(cnt).compute_grad_Laplacian('tskin');
        RHB_alongwind_segobj(cnt)= RHB_alongwind_segobj(cnt).compute_grad_Laplacian('u_algtraj_true','type','gradient');
        RHB_alongwind_segobj(cnt)= RHB_alongwind_segobj(cnt).compute_grad_Laplacian('u_algtraj_10N','type','gradient');
        

        
        
        
        
        % estimate the length scale of the SST front:
        
    end
   % close all
    %pause
end

%% get RHB crosswind segments: (wind align with SST gradient)
cnt = 0;
for iseg = 1:nseg
    if ~isempty(RHB_straight_segs(iseg).crx_wind)
        cnt = cnt+1;
        
        % 3. project wind on to the trajectory
        % initialize object:
        RHB_crxwind_segobj(cnt) = ATOMIC_dataProcess(RHB_straight_segs(iseg).crx_wind);
        RHB_crxwind_segobj(cnt).ATOMIC_platform = 'RHB';

        
        RHB_crxwind_segobj(cnt) = RHB_crxwind_segobj(cnt).project_wind_onto_trajectory('wspd_varn','wspd');   
        RHB_crxwind_segobj(cnt) = RHB_crxwind_segobj(cnt).project_wind_onto_trajectory('wspd_varn','wspd_10N');   
        
        
        % compute SST gradient and laplacian using (central difference):
        RHB_crxwind_segobj(cnt) = RHB_crxwind_segobj(cnt).compute_grad_Laplacian('tsea');
        RHB_crxwind_segobj(cnt) = RHB_crxwind_segobj(cnt).compute_grad_Laplacian('tskin');
        
        RHB_crxwind_segobj(cnt)= RHB_crxwind_segobj(cnt).compute_grad_Laplacian('u_algtraj_true','type','gradient');
        RHB_crxwind_segobj(cnt)= RHB_crxwind_segobj(cnt).compute_grad_Laplacian('u_algtraj_10N','type','gradient');
        

        % no need to orient data for the cross wind.

        % estimate the length scale of the SST front:
    
    end
    
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
svdataname = ['RHB_along_wind_segments_ready_for_wavelet_coherence_updated_' date '.mat'];   %_crit80m
save([svdatadir filesep svdataname],'RHB_alongwind_segobj','RHB_crxwind_segobj','RHB_straight_segs','RHB_seglens');



%load([svdatadir filesep svdataname]);

%% put data into the wavelet coherence toolbox.
% add wind_div with SST gradient, wind_div with SST laplacian  (wind div.
% is evaluated with the wind component projected along the RHB
% trajectories)
%varnlist = { 'wspd_10N', 'wdir','u_algtraj'};  % add the directly measured wind speed;

varnlist = {'u_algtraj_true', 'u_algtraj_10N'};   % use the true wind speed and 10-m neutral wind

% here, 'wspd' refer to the along-transect wind component. According to
% these matrix pair, feed in the relevant parameters and conduct the
% wavelet analysis. 
matrix_pairs = [{'wspd','SST'}; {'wspd_div', 'SST_grad'};  ...
               {'wspd_div','SST_laplacian'};{'wspd_div','dSST_dt'};];
           
% create different data and figure output directory and data file name;
svdata_subdir = [svdatadir filesep 'test_different_matrix_pairs'];
if ~exist(svdata_subdir, 'dir')
    mkdir(svdata_subdir);
end

% I need to specify the right label for each matrix pairs as well.  
SST_varn = 'tsea';
           
svfig = true;

for iv = 1 %:length(varnlist)
    wind_varn = varnlist{iv};
    
    loc_svdatadir = [svdata_subdir, filesep wind_varn];
    if ~exist(loc_svdatadir, 'dir')
        mkdir(loc_svdatadir);
    end
    
    loc_figsvdir = [loc_svdatadir filesep 'figs'];
    if ~exist(loc_figsvdir, 'dir')
        mkdir(loc_figsvdir);
    end
    
    % add data pair to test:
    wnddiv_varn = [wind_varn '_grad'];      % I need to compute this term
    SSTgrad_varn = [SST_varn '_grad'];
    SSTcurv_varn = [SST_varn '_curv'];
    dSSTdt_varn = [SST_varn '_tendency'];   % and this one;
    
    data_pair = [{wind_varn, SST_varn}; {wnddiv_varn , SSTgrad_varn}; ...
                 {wnddiv_varn, SSTcurv_varn};  {wnddiv_varn,dSSTdt_varn}; ];
    label_pair{1} = {'SST (^{\circ}C)', {'along-transect','wind magnitude (m/s)'}};
    label_pair{2} = {'SST gradient (^{\circ}C/km)', {'along-transect', 'wind divergence (s^{-1}x1000)'}};
    label_pair{3} = {'SST curvature (^{\circ}C/km^2)', {'along-transect', 'wind divergence (s^{-1}x1000)'}};
    label_pair{4} = { {'effective SST', 'rate-of-change (^{\circ}C/s x 1000)'}, {'along-transect', 'wind divergence (s^{-1}x1000)'}};
 
    
    for p = 1:size(data_pair,1)
        v1 = cell2mat(data_pair(p,2));    % SST related variables
        v2 = cell2mat(data_pair(p,1));    % wind related variables
        
        clear RHB_wtc_stat RHB_inCOI_stat RHB_outCOI_stat
        jj=0;
        for i = 1:length(RHB_alongwind_segobj)
            inobj = RHB_alongwind_segobj(i);     % how would the sampling resolution and remapping resolution change the results?
            
            if p==4 % compute the effective SST rate of change experienced by a parcle:
                % compute this variable on the fly for different
                % combination of SST and wind speed;
                inobj = inobj.compute_effective_SSTrate_of_change(SST_varn, wind_varn);
            end
                
            for iseg =1:length(inobj.Value)
                jj = jj+1;
                %validID = ~isnan(ts1);
                ts1 =inobj.Value(iseg).(v1);
                %ts1= ts1
                %ts2 = WG247_alongwind_segobj_remapped.Value(iseg).wind_speed;
                ts2 = inobj.Value(iseg).(v2);      % take the absolute values to remove direction?
                traj = inobj.Value(iseg).distance_axis - inobj.Value(iseg).distance_axis(1);
                time = inobj.Value(iseg).time;
                
                ts_pair(jj).ts1 = ts1;
                ts_pair(jj).ts2 = ts2;
                ts_pair(jj).traj = traj;
                ts_pair(jj).time = time;
                
                ATOMIC_platform = 'RHB';
                
                labelstr.ts1 = label_pair{p}{1};
                labelstr.ts2 = label_pair{p}{2}; %'along-traj.';
                
                %%% 0. establish the wavelet object for each glider:
                RHB_wtcobj = ATOMIC_Wavelet(ts1, ts2, traj, time, ATOMIC_platform, RHB_xres);
                
                if svfig
                    hfig0 = RHB_wtcobj.plot_data(labelstr);
                    xc_savefig(hfig0, loc_figsvdir, [ATOMIC_platform '_seg' num2str(i, '%2.2i') ...
                        'part-' num2str(iseg) '_' wind_varn '_'  v1 '_' v2 '.jpg'], [0, 0, 10 8]);
                end
                
                if p ==1
                    %%% 1. apply box-cox transform on the data to have normal distribution;
                    [RHB_wtcobj_trans, hfig1] = RHB_wtcobj.transform_to_normal_distribution(labelstr);
                    
                    if svfig
                        xc_savefig(hfig1, loc_figsvdir, [ATOMIC_platform '_seg' num2str(i, '%2.2i') ...
                            'part-' num2str(iseg) '_' wind_varn '_' v1 '_' v2 '_boxcox_transformed.jpg'], [0, 0, 10 8]);
                    end
                    
                    %
                    
                    %%% 2. apply the wave coherence toolbox:
                    
                    % how to make sure the x-axis is labelled correctly. 
                    [hfigs, RHB_wtc_stat(jj), RHB_outCOI_stat(jj), RHB_inCOI_stat(jj)] = RHB_wtcobj_trans.wavelet_coherence_toolbox;
                    
                else
                    % directly apply wavelet coherence:
                    [hfigs, RHB_wtc_stat(jj), RHB_outCOI_stat(jj), RHB_inCOI_stat(jj)] = RHB_wtcobj.wavelet_coherence_toolbox;
                end
                    
                
                    figure(hfigs.Number);
                    title([ATOMIC_platform ': ' datestr(time(1)) '~' datestr(time(end))]);    %UTC time;

                if svfig
                    xc_savefig(gcf, loc_figsvdir, [ATOMIC_platform '_seg' num2str(i, '%2.2i') ...
                        'part-' num2str(iseg) '_coherence_boxcox_transformed_' wind_varn '_' v1 '_' v2 '.jpg'], [0, 0, 10 8]);
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
        svdataname = ['RHB_wavecoherence_stat_' v1 '_and_' v2 '_xres' num2str(RHB_xres) 'km_' timestamp '.mat'];
        save([loc_svdatadir filesep svdataname],'RHB_wtc_stat','RHB_outCOI_stat','RHB_inCOI_stat','ts_pair');
        
    end
    
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
