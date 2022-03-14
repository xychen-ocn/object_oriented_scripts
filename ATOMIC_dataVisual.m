% ======================================================================= %
% class: ATOMIC_dataVisual
% purpose: visualize ATOMIC SST data:
% Author: Xuanyu Chen 
% Date: v0 @ 09/29/2021
% ======================================================================= %
classdef ATOMIC_dataVisual
    properties
        data struct
        colormap = 'parula'
        linewidth = 1.2
        fontsize = 14
        fignum = 1
        coastline_data='/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/satellite_SLA/NACoast_shoreline_h.mat';
        scasize = 20
        marker = 'o';
        linestyle = '-';
        color = 'k';
        markerEdgeColor ='none';
    end
    
    
    methods
       % -------------------- Function # 1: ------------------------- %
       % -----       initialize the object:             ----- %
       function obj = ATOMIC_dataVisual(indata, figppt)
           switch nargin
               case 1
                   obj.data = indata;
                   
               case 2
                   obj.data = indata;
                   ppts = fieldnames(figppt);
                   for i = 1:length(ppts)
                       obj.(ppts{i}) = figppt.(ppts{i});
                   end
                   
               case 0                   
                   data = struct();
                   % everthing else will be defalut.
           end
       end
           

       % --------------------- Function # 2: -------------------------- %
       % -----    map out the available measurements:             ----- %
       function [h, hb, ax] = collect_dataset_on_map(obj, cc_datan)
           % output figure handle and axis handle.
           h = figure(obj.fignum); hold on;
           % make use of coastline data:
           load(obj.coastline_data);
           plot(coast.lon, coast.lat, '-k', 'linewidth',1.5);
           patch(coast.lon, coast.lat, [0.5 0.5 0.5]);
           % plot the measurement track colorcoded by time of the
           % measurement.
           data = obj.data;
           cc_data = data.(cc_datan);
           if numel(cc_data)~=numel(data.lon)
               % AXBTs
               cc_data = mean(cc_data(1:10,:),1,'omitnan');
           end
           scatter(data.lon, data.lat, obj.scasize, cc_data, 'filled',...
               'marker', obj.marker,'markerEdgeColor',obj.markerEdgeColor);
           
           xlabel('Longitude (^{\circ}E)');
           ylabel('Latitude (^{\circ}N)');
           colormap(obj.colormap);
           hb = colorbar;
           %datetick(hb,'x','mmmdd HH','keepticks');
           %set(get(hb,'xlabel'),'String',)
           x0 = floor(min(data.lon)*1.01); xn = ceil(max(data.lon)*0.99);
           y0 = floor(min(data.lat)*0.99); yn = ceil(max(data.lat)*1.01);
           
           axis([x0, xn, y0, yn])
           set(gca,'fontsize',obj.fontsize);
           
            ax = gca;

           
           
       end
       

       % --------------------- Function # 2: -------------------------- %
       % -----    map out the available measurements:             ----- %
       function [hl] = add_trajectory(obj, lineppt)
           % output figure handle and axis handle.
           h = figure(obj.fignum); hold on;

           data = obj.data;
           hl = plot(data.lon(1:5:end), data.lat(1:5:end), 'linestyle',lineppt.linestyle, ...
               'linewidth',lineppt.linewidth, 'color',lineppt.color, 'marker', lineppt.marker);
           
           plot(data.lon(1), data.lat(1), 'p','color', lineppt.color, ...
               'markersize', lineppt.markersize, 'linewidth',lineppt.linewidth);
           
           plot(data.lon(end), data.lat(end), 'd','color', lineppt.color, ...
               'markersize', lineppt.markersize, 'linewidth',lineppt.linewidth);
           
           if lineppt.textflag
               text(mean(data.lon), mean(data.lat)*1.005,datestr(data.local_time(1),'mmmdd'),'fontsize',12,'fontweight','bold');
           end
           
           
           
           
       end

       
       % --------------------- Function # 2: -------------------------- %
       % -----    plot time series            ----- %
       function hl = plot_time_series(obj,sstname)
           
           h = figure(obj.fignum); hold on;
           data = obj.data;
           local_time = data.time - 4/24;
           sst = data.(sstname);
           hl = plot(local_time, sst);
           hl.LineStyle = obj.linestyle;
           hl.LineWidth = obj.linewidth;
           hl.Color = obj.color;
           
           datetick('x','HH:MM');
           ylabel('SST (^{\circ}C)');
           set(gca,'fontsize',obj.fontsize);
       end
       
       
       
       % --------------------- Function # 2: -------------------------- %
       % -----    plot time series            ----- %
       function hl = plot_variation_along_track(obj,sstname)
           
           h = figure(obj.fignum); hold on;
           data = obj.data;
           local_time = data.time - 4/24;
           sst = data.(sstname);
           hl = plot(local_time, sst);
           hl.LineStyle = obj.linestyle;
           hl.LineWidth = obj.linewidth;
           hl.Color = obj.color;
           
           datetick('x','HH:MM');
           ylabel('SST (^{\circ}C)');
           set(gca,'fontsize',obj.fontsize);
       end
       
       % --------------------- Function # 3: --------------------------- %
       % plot the along trajectory velocity:
       function hfig = plot_current_along_track(obj)
           hfig = figure(obj.fignum); hold on;
           data = obj.data;
           nseg = length(data)
           for i = 1:nseg
           traj = obj.data(i).distance_axis;
           subplot(2,nseg,i)%1+2*(i-1)
           yyaxis left
           plot(traj, data(i).tsea - mean(data(i).tsea,'omitnan'),'--r');hold on
           plot(traj, data(i).cspd,'-');
           hold on;
           yyaxis right
           plot(traj, data(i).cdir);
           
           subplot(2,nseg,i+nseg)
           cdir_cart = get_cartesian_direction(data(i).cdir,'Meteo');
           cu= data(i).cspd.*cosd(cdir_cart);
           cv = data(i).cspd.*sind(cdir_cart);
           quiver(data(i).lon, data(i).lat, fliplr(cu), fliplr(cv),'k');
           hold on;
           scatter(data(i).lon, data(i).lat, 30, fliplr(data(i).tsea),'filled');
           end
       end

       
        
    end
    
    
    
    
end