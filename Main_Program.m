% This script is used to map out the available SST sample from EUREC4A:
%
% data: 
%  - Ron Brown;
%  - SWIFTS
%  - Wave Glider:
%  - P3 IR radiometer:
%  - P3 AXBTs;


clear all; close all; clc;
%%%%%% 1. load in all the datasets %%%%%%
% - read each dataset from the source and separate the data by days. 
locdir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC';
rmtdir_psl = 'https://psl.noaa.gov/thredds/dodsC/Datasets/ATOMIC/data';
rmtdir = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/psl/atomic';
addpath([locdir filesep 'rhb/object_oriented_scripts']);
figsvdir = [locdir filesep 'rhb/object_oriented_scripts/Figs'];

datadir.RHB = [rmtdir_psl filesep 'rhb/nav_met_sea_flux'];
datadir.RHB = [locdir filesep 'rhb/data'];
datadir.SWIFTs = [locdir filesep 'SWIFTs_WG'];
datadir.WaveGlider = [locdir filesep 'SWIFTs_WG'];
datadir.P3_IR = [locdir filesep 'rhb/data'];
datadir.P3_AXBTs = datadir.P3_IR;

% datadir.P3_IR = [rmtdir filesep 'p3/remote_sensing'];
% datadir.P3_AXBTs = [rmtdir filesep 'p3/AXBT'];

% dataFN.RHB_all = 'rhb_daily_grouped_10min_data_0909latest.mat';
% dataFN.RHB_grouped = 'rhb_daily_grouped_10min_data_0909latest.mat';
dataFN.RHB = 'EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
dataFN.RHB = 'EUREC4A_ATOMIC_RHB_AllinOne_v1.3.mat';

% dataFN.SWIFTs = 'EUREC4A_ATOMIC_SWIFT*_All_v2.2.nc';
% dataFN.WaveGlider = 'EUREC4A_ATOMIC_WG*_All_v2.2.nc';
% dataFN.P3_IR = 'EUREC4A_ATOMIC_P3_Remote-sensing_*_v1.1.nc';
% dataFN.P3_AXBTs = 'EUREC4A_ATOMIC_P3_AXBT_v1.0.nc';

dataFN.P3_IR = 'EUREC4A_ATOMIC_P3_Flight-Level_*_v1.1.nc';

dataFN.SWIFTs = 'EUREC4A_ATOMIC_SWIFTs_AllinOne_v2.2.mat';
dataFN.WaveGlider = 'EUREC4A_ATOMIC_WaveGlider_AllinOne_v2.2.mat';
dataFN.P3_IR = 'EUREC4A_ATOMIC_P3_IR_AllinOne_v1.1.mat';
dataFN.P3_AXBTs = 'EUREC4A_ATOMIC_P3_AXBTs_AllinOne_v1.0.mat';


% the structure above will be used as to initialize the ATOMIC_dataIO
% class:

DataSet_Names = fieldnames(datadir);
NumDataSet = length(DataSet_Names);

for i = 1%:NumDataSet
    DSN = DataSet_Names{i};
    % prepare object properties to initialize objects:
    ppt.DataPath = datadir.(DSN);
    ppt.DataName = dataFN.(DSN);
    if contains(ppt.DataPath, 'xchen')
        ppt.DataLoc = 'local';
    else        
        ppt.DataLoc = 'opendap';
    end
    ppt.saveflag = true;
    ppt.DataTag = DSN;
    
    % initialize an ATOMIC data obj:
    data_obj = ATOMIC_dataIO(ppt);
    
    % use this object to either read a data from the server or load local
    % data:
    if contains(ppt.DataName, 'mat')
        % load in matlab data:
       data = data_obj.load_local_matdata;       
    else       
       % read netCDF first:
       data = data_obj.read_netCDF_to_matdata;     
    end
    
    eval([DSN '_ds = data.(DSN);']);
    clear data
    
    %pause
end

% group RHB data by (local time):
obj_tmp = ATOMIC_dataProcess(RHB_ds);
RHB_grped = obj_tmp.group_data_by_day_localtime;

% group SWIFTs data by (local time day)
SWIFTs_ds = SWIFTs_ds.SWIFTs;
floatIDs = fieldnames(SWIFTs_ds);
NumSWIFTs = length(floatIDs);

for i = 1:NumSWIFTs
    FN = floatIDs{i};
    obj_tmp = ATOMIC_dataProcess(SWIFTs_ds.(FN));    % temporary object;
    SWIFTs_grped.(FN) = obj_tmp.group_data_by_day_localtime;
end
% save this data:


% group Wave Glider by local day:
WaveGlider_ds = WaveGlider_ds.WaveGlider;
WGIDs = fieldnames(WaveGlider_ds);
NumWGs = length(WGIDs);

for i = 1:NumWGs
    FN = WGIDs{i};
    obj_tmp = ATOMIC_dataProcess(WaveGlider_ds.(FN));    % temporary object;
    WGs_grped.(FN) = obj_tmp.group_data_by_day_localtime;
end
% save this data

P3_AXBTs_ds = P3_AXBTs_ds.P3_AXBTs;
P3_AXBTs_ds = convert_Kelvin_to_Celcius(P3_AXBTs_ds, 'temperature');
obj_tmp = ATOMIC_dataProcess(P3_AXBTs_ds);
P3_AXBTs_grped = obj_tmp.group_data_by_day_localtime;


% convert P3 radiometer SST from K to °C:
P3_IR_ds = P3_IR_ds.P3_IR;
P3obj = ATOMIC_dataProcess(P3_IR_ds);
P3_IR_ds = convert_Kelvin_to_Celcius(P3obj,'SST_IR_est');




%%%%% 2. produce a map %%%%%%
% - i. showing all the data on all days
% - ii. showing the SST map on 1 day with all available measurements.
% - ?

% Note: Here I will need a class: ATOMIC_dataVisual (for showing SST map on
% each day. 

% days of interest: Jan 9 (WG launched); Jan 14(SWIFTs launched); Jan 22
% (SWIFTs recovered); Jan 30 (SWIFTs launched leg2); Feb 11 (end of leg of RHB);


%% on the same figure, add 6 SWIFTs traj:
clf;
%%%%% Wave Gliders %%%%%
figppt.fignum =2;
figppt.scasize = 30;
figppt.marker = 's';
figppt.markerEdgeColor = 'none';

figppt.colormap = 'jet';
cc_datan = 'sea_water_temperature';

trajppt.linewidth=1.4;
trajppt.linestyle='none';
trajppt.color=[0 0.4470 0.7410];
trajppt.markersize=10;
trajppt.marker='none';
trajppt.textflag  = false;

hl_WGs=[];
for ii = 1:NumWGs
    FN = WGIDs{ii};
    % if isfield(WGs_grped.(FN), DN)
    figobj = ATOMIC_dataVisual(WaveGlider_ds.(FN), figppt);
    [hfig, hbar, ax]=collect_dataset_on_map(figobj, cc_datan);
    hl_WGs(ii) = add_trajectory(figobj, trajppt);
    %  end
end
grid on
set(ax, 'XMinorGrid','on','XMinorTick','on');
set(ax, 'YMinorGrid','on','YMinorTick','on');
plot(WG247_alongwind_segobj

title('Wave Gliders');
set(get(hbar,'xlabel'),'String','SST (^{\circ}C)');
xc_savefig(gcf,'./','WGs_SST_traj_demo.jpg',[0 0 10 8 ]); 


%%%%% SWIFTs %%%%%
figppt.fignum= 1;
figppt.scasize = 25;
figppt.marker = '^';
figppt.markerEdgeColor = 'none';
figppt.colormap = 'jet';
cc_datan = 'sea_water_temperature';

trajppt.linewidth=1.2;
trajppt.linestyle='none';
trajppt.color=[0.5 0.5 0.5];
trajppt.markersize=10;
trajppt.marker='.';
trajppt.textflag  = false;

hl_SWF=[];
for ii = 1:NumSWIFTs
    FN = floatIDs{ii};
    %if isfield(SWIFTs_grped.(FN), DN)
    figobj = ATOMIC_dataVisual(SWIFTs_ds.(FN), figppt);
    [hfig, hbar, ax]=collect_dataset_on_map(figobj, cc_datan);
    hl_SWF(ii) = add_trajectory(figobj, trajppt);
    
    %end
end
 
    

%% For RHB, add the trajectories on the following days of interest.
DOIs ={'Jan08','Jan09','Jan10',...
       'Jan18',...
       'Jan22','Jan23','Jan24', ...       
       'Feb03','Feb04', ...
       'Feb06','Feb07','Feb08', ...
       'Feb10','Feb11','Feb12', ...
       'Jan14','Jan30'};
   %'Jan13','Jan14','Jan15',...   % jan14: deployment of swifts
   %'Jan29','Jan30', ...
   
% DOIs = {'Jan18','Feb04','Feb08'};
% DOIs = {'Jan23'};
%clf;
for i = 1:length(DOIs)
    DN = DOIs{i};
    
    %%%%% RHB %%%%%
    %clf;
    figppt.fignum = 2;
    figppt.scasize = 50;
    figppt.marker = 'o';
    figppt.colormap = 'jet';
    cc_datan = 'tskin';
    trajppt.linewidth=1.2;
    trajppt.linestyle='-';
    trajppt.color='k';
    trajppt.markersize=12;
    trajppt.marker='none';
    trajppt.textflag  = true;
    
 
   % figobj = ATOMIC_dataVisual(RHB_grped.(DN), figppt);
     figobj = ATOMIC_dataVisual(RHB_ds, figppt);
    [hfig, hbar, ax]=collect_dataset_on_map(figobj, cc_datan);
    %hl_RHB = add_trajectory(figobj, trajppt);
end
grid on
set(ax, 'XMinorGrid','on','XMinorTick','on');
set(ax, 'YMinorGrid','on','YMinorTick','on');
axis('equal');
xlim([-62, -50]);
ylim([12 16.5]);
title('RV Ronald H. Brown from Jan-09 to Feb-12 2020')
% title('RHB and SWIFTs');
% title('RHB, SWIFTs & WGs');
set(get(hbar,'xlabel'),'String','SST (^{\circ}C)');

xlim = get(ax,'xlim');
ylim = get(ax,'ylim');

set(ax,'xlim',xlim);
set(ax,'ylim',ylim);

xc_savefig(gcf,'./','RHB_SST_traj_demo.jpg',[0 0 10 8 ]); 

xc_savefig(gcf,'./','RHB_SWIFTs_SST_traj_demo.jpg',[0 0 10 8 ]); 
xc_savefig(gcf,'./','RHB_SWIFTs_WGs_SST_traj_demo.jpg',[0 0 10 8 ]); 
   
xc_savefig(gcf,'./','RHB_SWIFTs_SST_traj_demo_upstream_trade_leg1.jpg',[0 0 10 8 ]); 
xc_savefig(gcf,'./','RHB_SWIFTs_WGs_SST_traj_demo_upstream_trade_leg1.jpg',[0 0 10 8 ]); 



figppt.scasize = 30;
figppt.marker = 's';
figppt.markerEdgeColor = 'none';

figppt.colormap = 'jet';
cc_datan = 'sea_water_temperature';

trajppt.linewidth=1.4;
trajppt.linestyle='none';
trajppt.color=[0 0.4470 0.7410];
trajppt.markersize=10;
trajppt.marker='none';
trajppt.textflag  = false;

hl_WGs=[];
for ii = NumWGs
    FN = WGIDs{ii};
    % if isfield(WGs_grped.(FN), DN)
    figobj = ATOMIC_dataVisual(WaveGlider_ds.(FN), figppt);
    [hfig, hbar, ax]=collect_dataset_on_map(figobj, cc_datan);
    hl_WGs(ii) = add_trajectory(figobj, trajppt);
    %  end
end






    
    title([DN ' (local time)'],'fontsize',16);
    caxis([25.6, 27]);
    if ~isempty(hl_SWF) && ~isempty(hl_WGs)
        ids = find(hl_SWF>0);
        ids2 = find(hl_WGs>0);
        hline = [hl_RHB, hl_SWF(ids(1)), hl_WGs(ids2(1))];
        lgdstr = {'RHB','SWIFT','WG'};
    else
        if ~isempty(hl_SWF) && isempty(hl_WGs)
            ids = find(hl_SWF>0);
            hline = [hl_RHB, hl_SWF(ids(1))];
            lgdstr = {'RHB','SWIFT'};
        elseif isempty(hl_SWF) && ~isempty(hl_WGs)
            ids2 = find(hl_WGs>0);
            hline = [hl_RHB, hl_WGs(ids2(1))];
            lgdstr = {'RHB','WG'};
            
        else isempty(hl_SWF) && isempty(hl_WGs)
            hline = hl_RHB;
            lgdstr= 'RHB';
        end
    end
    
    
    
    %%%% P3 trajectory %%%%
    figppt.scasize = 20;
    figppt.marker = '*';
    trajppt.linewidth=1.1;
    trajppt.linestyle='--';
    trajppt.color=[0.9290 0.6940 0.1250];
    trajppt.markersize=10;
    
    if isfield(P3_IR_ds, DN)
        figobj = ATOMIC_dataVisual(P3_IR_ds.(DN), figppt);
        hl_P3 = add_trajectory(figobj, trajppt);
        hline = [hline hl_P3];
        lgdstr{end+1} = 'P3';
    end
    hlgd = legend(hline,lgdstr);

        
    set(get(hbar,'xlabel'), 'String','SST (^{\circ}C)')
    
    pause()
    
    
    figname = [DN '_SST_map_from_various_instruments_at_instrument_depthlevs.jpg']
    xc_savefig(hfig, figsvdir, figname, [0 0 10 8]);

    
end


DOIs ={'Jan09','Jan14','Jan22','Jan30','Feb11'};
% plot timeseries
DOIs =  {'Jan18','Jan23','Feb04','Feb08'};

for i = 1:length(DOIs)
    DN = DOIs{i};
    
    %%%%% RHB %%%%%
    clf;
    figppt.fignum = i;
    figppt.scasize = 50;
    figppt.marker = 'o';
    figppt.colormap = 'jet';
    cc_datan = 'tsea';
    figppt.linewidth=1.2;
    figppt.linestyle='-';
    figppt.color='k';
    
    
    figobj = ATOMIC_dataVisual(RHB_grped.(DN), figppt);
    hl_RHB = plot_time_series(figobj, cc_datan);
    
    figppt.color='b';
    figppt.linewidth=0.8;
    figppt.linestyle='-';
    figobj = ATOMIC_dataVisual(RHB_grped.(DN), figppt);
    plot_time_series(figobj, 'tskin');
    
    figppt.color='r';
    figobj = ATOMIC_dataVisual(RHB_grped.(DN), figppt);
    plot_time_series(figobj, 'tsea_ship');
    
    %      figppt.color='r';
    %      figobj = ATOMIC_dataVisual(RHB_grped.(DN), figppt);
    %      plot_time_series(figobj, 'tsea_prof_8m');
    legend({'tsea (0.05m)','tskin (interfacial)','tsea_{ship} (5m)'})
    figname = [DN '_SST_timeseries_from_compared_at_3depthlevs.jpg']
    xc_savefig(gcf, figsvdir, figname, [0 0 10 8]);
    
end
    
    %%%%% SWIFTs %%%%%
    figppt.scasize = 40;
    figppt.marker = '^';
    figppt.colormap = 'jet';
    cc_datan = 'sea_water_temperature';
    
    figppt.linewidth=1.2;
    figppt.linestyle='--';
    colors = parula(12);
    
    hl_SWF=[];
    for ii = 1:NumSWIFTs
        FN = floatIDs{ii};
        if isfield(SWIFTs_grped.(FN), DN)
            figppt.color=colors(ii,:);

            figobj = ATOMIC_dataVisual(SWIFTs_grped.(FN).(DN), figppt);
            hl_SWF(ii) = plot_time_series(figobj, cc_datan);

        end
    end
 
    

    %%%%% Wave Gliders %%%%%
    figppt.scasize = 40;
    figppt.marker = 's';
    figppt.colormap = 'jet';
    cc_datan = 'sea_water_temperature';
    
    figppt.linewidth=1.4;
    figppt.linestyle='-.';
    colors = spring(5);

    hl_WGs=[];
    for ii = 1:NumWGs
        FN = WGIDs{ii};
        if isfield(WGs_grped.(FN), DN)
            figppt.color = colors(ii*2,:);
            figobj = ATOMIC_dataVisual(WGs_grped.(FN).(DN), figppt);
             hl_WGs(ii) = plot_time_series(figobj, cc_datan);
        end
    end

    
    
    title([DN ' (local time)'],'fontsize',16);
    caxis([25.6, 27]);
    if ~isempty(hl_SWF) && ~isempty(hl_WGs)
        ids = find(hl_SWF>0);
        ids2 = find(hl_WGs>0);
        hline = [hl_RHB, hl_SWF(ids(1)), hl_WGs(ids2(1))];
        lgdstr = {'RHB','SWIFT','WG'};
    else
        if ~isempty(hl_SWF) && isempty(hl_WGs)
            ids = find(hl_SWF>0);
            hline = [hl_RHB, hl_SWF(ids(1))];
            lgdstr = {'RHB','SWIFT'};
        elseif isempty(hl_SWF) && ~isempty(hl_WGs)
            ids2 = find(hl_WGs>0);
            hline = [hl_RHB, hl_WGs(ids2(1))];
            lgdstr = {'RHB','WG'};
            
        else isempty(hl_SWF) && isempty(hl_WGs)
            hline = hl_RHB;
            lgdstr= 'RHB';
        end
    end
    hlgd = legend(hline,lgdstr);

        
    
    pause()
    
    
    figname = [DN '_SST_timeseries_from_various_at_instrument_depthlevs.jpg']
    xc_savefig(gcf, figsvdir, figname, [0 0 10 8]);

    
end


% is SST gradient in skin-temperature different from that in the sea-snake
% temperature? (no)
%%%%% RHB %%%%%
%clf;
figppt.fignum = i;
figppt.scasize = 50;
figppt.marker = 'o';
figppt.colormap = 'jet';
cc_datan = 'tsea';
figppt.linewidth=1.2;
figppt.linestyle='-';
figppt.color='k';


figobj1 = ATOMIC_dataVisual(RHB_grped.(DN), figppt);
hl_RHB = plot_time_series(figobj1, 'tsea');

figppt.color = 'r';
figobj2 = ATOMIC_dataVisual(RHB_grped.(DN), figppt);
hl_RHB2 = plot_time_series(figobj2, 'tsea_prof_8m');



swift_obj = ATOMIC_dataVisual(SWIFTs_ds.SWIFT16, figppt);
plot_time_series(swift_obj, 'sea_water_temperature');

colors = spring(5);
figppt.fignum =10;
figppt.color=colors(2,:);
WG_obj = ATOMIC_dataVisual(WaveGlider_ds.WG245, figppt);
hl_wg245=plot_time_series(WG_obj, 'sea_water_temperature');
hl_wg245.Color = colors(2,:);

figppt.color=colors(4,:);
WG_obj = ATOMIC_dataVisual(WaveGlider_ds.WG247, figppt);
hl_wg247=plot_time_series(WG_obj, 'sea_water_temperature');
datetick('x','mmdd');
legend([hl_wg245, hl_wg247],{'WG245','WG247'})

figppt.scasize = 20;
figppt.marker = 's';
figppt.colormap = 'jet';
cc_datan = 'SST_IR_est';
figobj = ATOMIC_dataVisual(P3_IR_ds.(DN), figppt);
[hfig, hbar, ax]=collect_dataset_on_map(figobj, cc_datan);


figppt.scasize = 40;
figppt.marker = '<';
figppt.colormap = 'jet';
cc_datan = 'temperature';
figobj = ATOMIC_dataVisual(P3_AXBTs_grped.(DN), figppt);
[hfig, hbar, ax]=collect_dataset_on_map(figobj, cc_datan);


title([DN ' (local time)'],'fontsize',16);
caxis([25.6, 26.6]);






%% Part 3:
%%%%% 3. find SST gradient at mesoscale and save the segments found on each day   %%%%%
% input: variable name for SST. 
% output: a matlab structure of the segment.
% 
% structure: 
% -- loop 1: each datset;
% -- loop 2: each day (either local time or UTC) or?
% -- 
SWIFT_sstname = 'sea_water_temperature';
RHB_sstname = 'tskin';
spatial_res = 2;                % units km 
wvlen_cutoff = 20;              % units km 
figInfo.checkflag= true;
figInfo.saveflag = true;
figInfo.fignum = 10;


%% to be tested: (TO-DO!!!)
for i = 1:NumSWIFTs
    SWFN = floatIDs{i};
    % establish the object:
    SWFobj = ATOMIC_dataProcess(SWIFTs_ds.(SWFN));
    SWIFTs_segs.(SWFN) = SWFobj.find_mesoscale_SSTgrad(SWIFT_sstname, ...
        spatial_res, wvlen_cutoff, figInfo);
end

for i =1:NumWGs
    WGFN = WGIDs{i};
    WGobj = ATOMIC_dataProcess(WaveGlider_ds.(WGFN));
    WaveGlider_segs.(WGFN) = WGobj.find_mesoscale_SSTgrad(SWIFT_sstname, ...
        spatial_res, wvlen_cutoff, figInfo);
end

% cheated here: mannually load the rhbdates in:
DOIs = rhbdates.strong_SSTvar;
for id = 1 %:length(DOIs) .   pass test for the first front. (need to decide whether or not to apply low pass filter on the wind).
    DN = datestr(DOIs(id), 'mmmdd');
    RHB_obj = ATOMIC_dataProcess(RHB_grped.(DN));
    RHB_obj.ATOMIC_platform = 'RHB';
    %RHB_obj.dates = DN;
    [RHB_segs, RHB_obj] = RHB_obj.find_mesoscale_SSTgrad(RHB_sstname, ...
        spatial_res, wvlen_cutoff, figInfo);
end


    

%%%% 4. generate lagrangian trajectory for each SST gradient segment %%%%%
% 1. collect all the segments:
% 2. generate input to compute lagrangian trajectory based on the segment
% locations

% SWIFTs_segs.generate_lagrangian_trajectory()
% WaveGlider_segs.generate_lagrangian_trajectory()
% RHB_segs.generate_lagrangian_trajectory()
RHB_segs_obj = ATOMIC_LagrangianTraj();
RHB_segs_obj.ATOMIC_platform = 'RHB';
RHB_segs_obj.TrajModel = 'LagTraj';
RHB_segs_obj.SimDuration= 6;
RHB_segs_obj.TrajLevs= 700;

RHB_seg_fields = fieldnames(RHB_segs);
num_RHBdays = length(RHB_seg_fields);
for id = 1:num_RHBdays
    DN = RHB_seg_fields{id};
    RHB_segs_obj.data = RHB_segs.(DN);
    RHB_segs_obj.generate_input_for_trajectory_simulation
end

% this will need to be execute from the linux machine.. (or, I download
% data locally..)
RHB_segs_obj.SimOutputDir = RHB_segs_obj.SimInputDir;
RHB_lagtraj = RHB_segs_obj.get_output_from_trajectory_simulation();

%% similarily:
SWIFTs_segs.generate_input_for_trajectory_simulation()
WaveGlider_segs.generate_input_for_trajectory_simulation()

SWIFTs_lagtraj = get_output_from_trajectory_simulation();
WaveGliders_lagtraj = get_output_from_trajectory_simulation();


%%
%%%% 5. construct cloud scene with its center moving along the lagrangian
%%%% trajectory, then compute various statistics                     %%%%%%
% this can be coded as well. 

%%%%

% add one more step here: prepare input for extracting the GOES16 data
SWIFTs_GOES16_input = SWIFTs_seg.prep_SateData_Retrieval; % not sure what I was trying to do here...
% multiple segments..

box_width = 10;        % units: degree (box width will be based on the trajectory;
time_window = 6;       % units: hr (length of LagTraj / average boundary layer wind speed. 
time_window = round(box_width*111 / (10x3.6));
% these two values can be computed within the function as well.

RHB_GOES16_obj = ATOMIC_GOES16;
RHB_GOES16_obj.data = datasub;    % manually load this from the right data
RHB_GOES16_obj.SateSource = '2km_10min_fulldisk';
RHB_GOES16_obj.SateDir = '';
RHB_GOES16_obj.SateChannel = 'VIS';
RHB_GOES16_obj.type = 'SateData';
RHB_GOES16_obj.DateLST = 'Jan09';


RHB_GOES_IR = RHB_GOES16_input.extract_GOES16_data('IR', boxwidth, time_window);
RHB_GOES_VIS = RHB_GOES16_input.extract_GOES16_data('VIS', boxwidth, time_window);




region_mask = SWIFTs_segs.construct_examination_region();

cloud_metrics_IR = SWIFTs_GOES_IR.compute_cloud_metrics(region_mask);
cloud_metrics_VIS = SWIFTs_GOES_VIS.compute_cloud_metrics(region_mask);

cloud_metrics_IR.visualize();
cloud_metrics_VIS.visualize();








