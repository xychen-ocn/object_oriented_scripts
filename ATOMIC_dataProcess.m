% ======================================================================= %
% class: ATOMIC_dataProcess
% purpose: process ATOMIC data
% Author: Xuanyu Chen 
% Date: v0 @ 09/29/2021
% ======================================================================= %

classdef ATOMIC_dataProcess
    properties
        Value struct
        FigOutput_rootdir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/object_oriented_scripts/Figs';
        ATOMIC_platform
    end
    
    methods
        %%%% ===================  function #0  ===================== %%%%
        %%%% initialize the object:
        function obj = ATOMIC_dataProcess(inval)
            if nargin > 0 
                if isstruct(inval)
                    obj.Value = inval;
                end
            else
                obj.Value = struct();  % initialize with an empty structure.
            end
        end
        
        
        %%%% ==================   function #1  ====================== %%%%        
        %%%%% : regroup data by day (local time or UTC)
        function ds_grouped = group_data_by_day_localtime(obj)
            
            ds = obj.Value;
            tmin = floor(min(ds.time))-1;
            tmax = ceil(max(ds.time))+1;
            tarry = tmin:1:tmax;
            
            ds.local_time = ds.time - 4/24; 
            [NY,NX] = size(ds.local_time);
            if NY~=1
                NT = NY;
            else
                NT = NX;
            end
            [Y, tedges] = discretize(ds.local_time, tarry);
            
            ds_fieldn = fieldnames(ds);
            
            % now group data according to dates:
            %obsv_timenum = [];
            cnt =0;
            for i = 1:length(tedges)
                ids = find(Y==i);
                DN = datestr(tedges(i),'mmmdd');
                if ~isempty(ids)
                    cnt = cnt+1;
                    for iv = 1:length(ds_fieldn)
                        varn = ds_fieldn{iv};
                        % need to do this according to the dimension of the
                        % variables:
                        [NY, NX] = size(ds.(varn));
                        if NY==NT
                            if NX==1
                                ds_grouped.(DN).(varn) = ds.(varn)(ids);
                            else
                                ds_grouped.(DN).(varn) = ds.(varn)(ids,:);
                            end
                        elseif NX==NT
                            if NY ==1
                                ds_grouped.(DN).(varn) = ds.(varn)(ids);
                            else
                                ds_grouped.(DN).(varn) = ds.(varn)(:,ids);
                            end
%                             ds_grouped.(varn) = ds.(varn);
                        end
                        
                    end
                    %obsv_timenum(cnt) = tedges(i);
                end
            end

            
        end
        
        
         %%%% ==================   function #1  ====================== %%%%        
        %%%%% : regroup data by day (local time or UTC)
        function ds_grouped = group_data_by_inputmask(obj, mask)
            
            ds = obj.Value;
           
            ds.local_time = ds.time - 4/24; 
            [NY,NX] = size(ds.local_time);
            if NY~=1
                NT = NY;
            else
                NT = NX;
            end
        
            ds_fieldn = fieldnames(ds);
            
            % now group data according to dates:
            %obsv_timenum = [];
            ids = mask;
            for iv = 1:length(ds_fieldn)
                varn = ds_fieldn{iv};
                % need to do this according to the dimension of the
                % variables:
                [NY, NX] = size(ds.(varn));
                if NY==NT
                    if NX==1
                        ds_grouped.(varn) = ds.(varn)(ids);
                    else
                        ds_grouped.(varn) = ds.(varn)(ids,:);
                    end
                elseif NX==NT
                    if NY ==1
                        ds_grouped.(varn) = ds.(varn)(ids);
                    else
                        ds_grouped.(varn) = ds.(varn)(:,ids);
                    end
                    %                             ds_grouped.(varn) = ds.(varn);
                end
                
            end
            %obsv_timenum(cnt) = tedges(i);
               
%             obj_new = obj;
%             obj_new.Value = ds_grouped;

            
        end
        
        %%%% ================= function ===================== %%%%
        %%%% purpose: make sure data doesn't have NaN values in the SST,
        %%%% location and .
        %%%% 
        function obj_new = drop_bad_records(obj)
            % this function is mainly for wave glider 245 that contains
            % NaN;
            data = obj.Value;
            
            valid_idx = find(isnan(data.lon)==0);
            edgeIDs = find(diff(valid_idx)~=1); % the separation loc;
            if length(valid_idx)<length(data.lon)
               
                disp('droping bad (NaN) records');
                if ~isempty(edgeIDs)
                dim1Dist = [edgeIDs(1); diff(edgeIDs);length(valid_idx)- edgeIDs(end)];
                else 
                    dim1Dist = [length(valid_idx)];
                end
                
                fieldn = fieldnames(data);
                for i = 1:length(fieldn)
                    FN = fieldn{i};
                    val = data.(FN);
                    % the edgeIDs store the cutting point of the data record
                    if iscolumn(data.lon)
                        tmp= mat2cell(val(valid_idx),dim1Dist, 1);
                    else
                        tmp = mat2cell(val(valid_idx), 1, dim1Dist);
                    end
                    
                    
                    % store data in different ways:
                    for j =1:length(tmp)
                        data_new(j).(FN) = tmp{j};
                    end
                end
                
                
                    
                
                obj_new = obj;
                obj_new.Value = data_new;
                
            else
                obj_new = obj;
            end
            
        end
            
        %%%% =================   function #2   ===================== %%%%
        %%%% Purpose: compute distance traveled by the ship/instrument
        function obj = compute_distance_travelled(obj)
            % input is expected to be free of NaN..
            data = obj.Value;     % this is a matlab structure;
            
            nseg = length(data);
            for iseg=1:nseg
                
                dlat = diff(data(iseg).lat);
                dlon = diff(data(iseg).lon);
                
                % average latitude between to consecutive ship locations:
                mid_lat = 0.5*(data(iseg).lat(1:end-1) + data(iseg).lat(2:end));
                xdist = abs(dlon) .* 111.*cosd(mid_lat);  % units: km
                ydist = abs(dlat) .* 111;
                
                % construct trajectory coordinate, (0,0) at the first location of the
                % ship.
                %
                if isrow(xdist)
                    traj_x = [0, cumsum(xdist)];
                    traj_y = [0, cumsum(ydist)];
                else
                    traj_x = [0; cumsum(xdist)];
                    traj_y = [0; cumsum(ydist)];
                end
                traj = sqrt(traj_x.^2 + traj_y.^2);
                
                data(iseg).traj = traj;
                
                % add computation of the course of direction:
                dy = dlat.* 111E3;
                dx = dlon.* 111E3.*cosd(mid_lat);
                
                trajdir_cart_mid = atan2(dy,dx).*180/pi;
                
                traj_mid = 0.5*(traj(1:end-1)+traj(2:end));
                
                nanmask = isnan(trajdir_cart_mid);
                trajdir_cart_mid(nanmask)=0;
                
                data(iseg).trajdir_cart = interp1(traj_mid, trajdir_cart_mid, traj, 'linear','extrap');
                
                
            end
            
            % extrapolate the direction for the last data point:
            
            
            obj.Value = data;
        end
        
        
        %%%% =================   function #2b   ===================== %%%%
        %%%% Purpose: compute distance traveled by the ship/instrument
        function obj = compute_ship_U10_misalignment(obj)
            
            data = obj.Value;     % this is a matlab structure;
            
            dlat = diff(data.lat);
            dlon = diff(data.lon);
            
            % find ship trajectory direction:
            mid_lat = 0.5*(data.lat(1:end-1) + data.lat(2:end));
            xdist = abs(dlon) .* 111.*cosd(mid_lat);  % units: km
            ydist = abs(dlat) .* 111;
            
            shipdir_cartesian = atan2(ydist, xdist)*180/pi;         % units: degree
            shipdir_cartesian(shipdir_cartesian<0) = shipdir_cartesian(shipdir_cartesian<0) + 360;  % make range to (0, 360)
            
            xmid = 0.5*( data.traj(1:end-1) + data.traj(2:end) ); 
            
            % convert the wind direction from the meteorolgoical convention
            % to cartesian coordinate:
            wind_direction_cartesian = get_cartesian_direction(data.wdir, 'Meteo');
            winddir_cartesian_midloc = interp1(data.traj, wind_direction_cartesian, xmid);
            
            ship_wind_angle = abs(shipdir_cartesian - winddir_cartesian_midloc);
            
            data.ship_wind_angle = ship_wind_angle;
            data.xmid = xmid;
            data.winddir_cart_xmid = winddir_cartesian_midloc;
            data.shipdir_cart = shipdir_cartesian;
            
            obj.Value = data;
            
            
          
        end
        
        
        %%%% =================   function #2b   ===================== %%%%
        %%%% Purpose: project data onto the RHB or wave glider trajectory. 
        function  [u_algtraj, v_crosstraj] = project_wind_on_trajectory(obj, prepped_data)
            
            data = obj.Value;
            
            % put x axis along trajectory and downwind.
            ang_diff = data.shipdir_cart - data.winddir_cart_xmid;
            against_locs = abs(ang_diff)>90;       % ship against wind:
            xaxis_dir = data.shipdir_cart;
            xaxis_dir(against_locs) = xaxis_dir(against_locs)+180;
            xaxis_dir(xaxis_dir>360) = xaxis_dir(xaxis_dir>360) - 360;
            
            ang_mid = xaxis_dir - data.winddir_cart_xmid;                  % this angle is 1 element less than the usual, at mid point of the two sample points.
            ang = interp1(data.xmid, ang_mid, data.traj,'linear','extrap');
            mask = prepped_data.moving_flag;
            ang_interp = interp1(data.traj(mask), ang(mask), prepped_data.dist_equal);
            
%             u_algtraj = data.U10N_interp .* cosd(ang_interp);                        % will always >0 , in the downstream direction
%             v_crosstraj = data.U10N_interp .* sind(ang_interp);                      % will change sign, southeasterly wind: (>0), northeasterly wind (<0);
            
            % low_passed: 
            u_algtraj = data.U10N_lowpassed .* cosd(ang_interp);                        % will always >0 , in the downstream direction
            v_crosstraj = data.U10N_lowpassed .* sind(ang_interp);                      % will change sign, southeasterly wind: (>0), northeasterly wind (<0);
           
            
        end
        
        
        %%%% =================   function #3    ==================== %%%%            
        %%%% Purpose: compute SST along-track gradient
        function SST_grad = compute_SST_gradient(obj)
            % refer to Ullman, D. S., & Cornillon, P. C. (2000). Evaluation of Front Detection Methods for Satellite-Derived SST Data Using In Situ Observations, Journal of Atmospheric and Oceanic Technology, 17(12), 1667-1675.
            % empirical though. 
            
            
            
            
        end
        
        
        %%%% =================   function #3    ==================== %%%%
        function moving_flag = get_ship_speed_mask(obj, spd_thres)
            data = obj.Value;
            if strcmpi(obj.ATOMIC_platform, 'RHB')
                moving_flag = data.sog>spd_thres;
            else
                moving_flag = true(size(data.time));
            end
        end
        
        
        
        
        %%%% ==================   function #4      ================= %%%%        
        %  find SST gradient at mesoscale
        function [seg, obj] = find_mesoscale_SSTgrad(obj, sstname, spatial_res, wvlen_cutoff, figInfo)
            
            % compute trajectory:
            obj = compute_distance_travelled(obj);
            
            % compute ship-wind misaglinment angle to determine whether the
            % ship is moving against the wind or not.
            obj = compute_ship_U10_misalignment(obj);

            
            % prepare the input SST and related quality to be on a
            % equally-space distance axis.
            % I might need to add a function that translate data into a
            % equally spaced axis first.
            [prepped_data, hfig_check, obj] = obj.prepare( sstname, spatial_res, wvlen_cutoff, figInfo);
            
            % The following function is a standalone function; (can be
            % improved to reduce the number of input values.   
            dataIn = obj.Value;
            seg = find_fronts_and_compute_gradients_OOP(dataIn, prepped_data,  wvlen_cutoff); 
            
            
            [u_algtraj, v_crosstraj] = project_wind_on_trajectory(obj,prepped_data);
            % compute the along track divergence and curl : du/dx
            %mask = prepped_data.moving_flag;
            
            % compute in the direction of the wind...
            wind_div_alongtraj = gradient(u_algtraj, fliplr(prepped_data.dist_equal));    % units: m/s per km . (using central difference)
            
            % calculate alongtrack wind curl: dv/dx
            wind_curl_alongtraj = gradient(v_crosstraj, fliplr(prepped_data.dist_equal));
            
            
            dataIn.wind_div_alongtraj_moving = wind_div_alongtraj;
            dataIn.wind_curl_alongtraj_moving = wind_curl_alongtraj;
            
            obj.Value = dataIn;
            
            % make plot for sanity check:
            if figInfo.checkflag
                sanity_check_results;
            end
            
            
            %% TO-DO!
            %----   Nested function here for plotting results out for
            %----      examination ----%% 
            % plot time series out instead (time-series of SST, wind speed,
            % and wind direction)
            function sanity_check_results
                fignum = figInfo.Number;
                figure(fignum); hold on;
                % plot both the segment identified and other properties
                % associated with the segments:
                
                moving_flag = prepped_data.moving_flag;
                
                subplot(2,2,[4]);
                hold on;
                plot(dataIn.lon, dataIn.lat,'--b');
                plot(dataIn.lon(1), dataIn.lat(1),'*k','markersize',13);
                plot(dataIn.lon(end), dataIn.lat(end),'dk','markersize',13);
                
                % overlay the surface wind direction:
                wndir = get_cartesian_direction(dataIn.wdir, 'Meteo');
                u = dataIn.wspd_10N(moving_flag).* cosd(wndir(moving_flag));
                v = dataIn.wspd_10N(moving_flag).* sind(wndir(moving_flag));
                
                if mode(wndir)>180
                    scale =-1;
                else
                    scale=1;
                end
                scatter(dataIn.lon(moving_flag), dataIn.lat(moving_flag)+ scale*0.025, 20, dataIn.wspd_10N(moving_flag),'filled');
                hb2=colorbar;
                set(get(hb2,'xlabel'),'string','m/s');
                %datetick(hb2,'x','mmmdd-hhZ','keepticks');
                caxis([median(dataIn.wspd_10N)-2*std(dataIn.wspd_10N), ...
                    median(dataIn.wspd_10N)+2*std(dataIn.wspd_10N)])
                %     quiver(dataIn.lon(1:3:end), dataIn.lat(1:3:end)-0.025, ...
                %         u(1:3:end), v(1:3:end), 0.5,'k');
                rhb_lon = dataIn.lon(moving_flag);
                rhb_lat = dataIn.lat(moving_flag);
                
                quiver(rhb_lon(1:3:end), rhb_lat(1:3:end)+ scale*0.025, ...
                    u(1:3:end), v(1:3:end), 0.5,'k');
                axis('square')
                xlabel('longitude');
                ylabel('latitude');
                set(gca,'fontsize',14)
                
                
                
                % make a along track plot for both the SST and wind speed
                % (& direction)
                % 1. plot SST detrend (change in SST)
                % 2. plot wind speed (change)
                % 3. compute along track convergence, and curl
                SST_detrend = prepped_data.SST_lowpassed_dtr;
                subplot(2,2,2);
                yyaxis left
                hl(1)=plot(prepped_data.dist_equal, fliplr(SST_detrend),'-b','linewidth',1.2);
                ylabel('SST (^{\circ}C)');
                set(gca,'fontsize',14);
                hold on;
                %hl(3)=plot(prepped_data.dist_equal, dataIn.wind_div_alongtraj_moving,'--k', 'linewidth',1.2);
                %
                %                 % calculate alongtrack wind curl:
               % hl(4)=plot(prepped_data.dist_equal, dataIn.wind_curl_alongtraj_moving,'-.k','linewidth',1.2);

                ylim([-0.5, 0.5])
                set(gca,'ytick',[-0.5:0.1:0.5]);
                
                yyaxis right
                wspd_demean = dataIn.U10N_interp - mean(dataIn.U10N_interp,'omitnan');
                wspd_detrend = detrend(dataIn.U10N_interp, 1);
                hl(2) = plot(prepped_data.dist_equal, fliplr(dataIn.U10N_interp), '-','linewidth',1.2); hold on
                ylabel('U_{10N} (m/s)');
                xlabel('distance traveled (km)');
                set(gca,'fontsize',14);
                
                % add alongtrack wind convergence:
                 %traj_moving = dataIn.traj(moving_flag);
                grid on
                
%                 
                hlgd = legend(hl,{'SST detrended','U_{10N}'}); %,'along traj div.'
                set(hlgd, 'loc','southwest');
                % save figure;
                figsvdir = [obj.FigOutput_rootdir filesep obj.ATOMIC_platform filesep 'find_SSTgrad'];
                if ~exist(figsvdir, 'dir')
                    mkdir(figsvdir);
                end
                DN = datestr(dataIn.local_time(1),'mmmdd');
                if figInfo.save
                figname = [DN '_SSTgrad_segment_v0_lowpassed.jpg']
                xc_savefig(gcf, figsvdir, figname, [0 0 14 10]);
                end
            end
            
        end
        
       
        %%% =================== function #6 ====================== %%%
        function obj_new = map_to_distance_axis(obj, spatial_res)
            data = obj.Value;
            
            for iseg = 1:length(data)
                if strcmp(obj.ATOMIC_platform , 'RHB')
                    spd_thres = 1;      % 1m/s-> 3.6km/hr;
                    moving_flag = get_ship_speed_mask(obj, spd_thres);
                else
                    moving_flag = true(size(data(iseg).traj));
                end
                traj = data(iseg).traj(moving_flag);
                dist_equal = [ceil(min(traj)):spatial_res:floor(max(traj))];
                
                % interpolation for variaous values:
                fieldn = fieldnames(data(iseg));
                for i = 1:length(fieldn)
                    FN = fieldn{i};
                    dimcrit = size(data(iseg).(FN)) == size(data(iseg).lon);
                    
                    if dimcrit
                        tmp_val = data(iseg).(FN)(moving_flag);
                        data_interp(iseg).(FN) = interp1(traj, tmp_val, dist_equal);
                    end
                    % skip variables that do not share the same dimension as
                    % lon.
                end
                data_interp(iseg).distance_axis = dist_equal;
                
                
                
                
                
            end
            obj_new = obj;
            obj_new.Value = data_interp;
            
        end
        
        
        %%%% ===================   function #5     =================== %%%%
        function [prepped_data, hfig, obj] = prepare(obj, sstname, spatial_res, wvlen_cutoff, ...
                                   figInfo)
        % Purpose: prepare SST record to get ready for front identification;
        % Steps: 1. interpolate SST record to an equally spaced record
        %        2. low-pass filter the SST record
        %        3. find peaks from the SST record
        %        4. 

        data = obj.Value;
        
        DN = datestr(data.local_time(1),'mmmdd');
        % the following code will need to be modified:, but it is a start.
         %%%%%%%%--> step 1:
         spd_thres = 1;      % 1m/s-> 3.6km/hr;
         moving_flag = get_ship_speed_mask(obj, spd_thres);
         sst = data.(sstname)(moving_flag);
         
         % map data to the distance axis:
         traj = data.traj(moving_flag);
         dist_equal = [0:spatial_res:round(max(traj))];                         % units: km
         SST_interp = interp1(traj, sst ,dist_equal,'linear', 'extrap');
         U10_interp = interp1(traj, data.wspd_10N(moving_flag),  dist_equal, 'linear','extrap');
         
         % the wind direction is likely not sensitive to the speed of the
         % ship(?)
         %data.ship_wind_angle_interp  = interp1(data.xmid, data.ship_wind_angle, dist_equal);
         
         
         %%%%%%%%--> step 2:
         %%% -- detrend:
         SST_interp_dtr = detrend(SST_interp,1);                                % detrend (I don't need this.)
         trend_linear = SST_interp - SST_interp_dtr;
         
         %U10N_interp_dtr = detrend(U10_interp,1);
         
      
         %%% -- low-pass filter:
         % wpass: normalized passband frequency: pi rad/sample.
         % filter out wavelength smaller than 20km;
         % wavenumber: 2*pi rad/20km;  resolution: 2km per sample
         % rad per sample = wavenumber_cutoff * resolution
         wn_cutoff = 2/wvlen_cutoff ;  % pi rad/km
         wpass = wn_cutoff * spatial_res;
         SST_lowpassed_dtr = lowpass(SST_interp_dtr, wpass);   % retain frequency lower than the cutoff;
         SST_lowpassed = SST_lowpassed_dtr+trend_linear;
         
         U10_lowpassed = lowpass_filter(U10_interp, wvlen_cutoff, spatial_res);
         
         %%%%%%%%--> step 3:
         
         minpeakdist = wvlen_cutoff/spatial_res;
         
         [pks.val, pks.loc, pks.w, pks.p] = findpeaks(SST_lowpassed, dist_equal, 'MinPeakDistance',minpeakdist, 'MinPeakProminence',0.05);
         [trghs.val, trghs.loc, trghs.w, trghs.p] = findpeaks(-(SST_lowpassed), dist_equal, 'MinPeakDistance',minpeakdist,'MinPeakProminence',0.01);
         
         %%%%%%% --> Step 4: continuous wavelet transform:
         [wt,fn]=cwt(SST_interp,'morse');
         
         % -- produce a sanity check figure:
         mag_wt = abs(wt);
         wvlen = 1./fn .* spatial_res;    % 1cycle/ (cycle per sample) * spatial resolution (km/sample) = km
         
         [DD, WL] = meshgrid(dist_equal, wvlen);
         BW = imregionalmax(abs(wt));
         % take the maximum that has a wavelength > 50km
         C1 = WL>25;
         C2 = mag_wt>0.75*max(abs(wt(:)));     % magnitude criterion set to 75% of the maximum magnitude.
         cond = (C1 & C2 & BW);
         
         % lengthscale of the local maximum :
         wvlen_locmax = WL(cond);
         xloc_locmax = DD(cond);
         
         %% plot
         if figInfo.checkflag
             hfig = figure(figInfo.fignum);clf;
             %%========== signal ========= %
             subplot(2,2,1)
             %yyaxis left
             plot(dist_equal,SST_interp,'linewidth',1.1)
             hold on
             plot(traj, sst,'.b');
             plot(dist_equal,trend_linear,'--k');
             xlabel('km');
             ylabel('Celcius');
             set(gca,'fontsize',14);
             
             %yyaxis right
             %plot(dist_equal,detrend(SST_interp,0),'o');
             hold on
             plot(pks.loc, pks.val,'*m','linewidth',1.2, 'markersize',12);
             plot(trghs.loc, -trghs.val,'*b','linewidth',1.2, 'markersize',12);
             plot(dist_equal, SST_lowpassed,'-r','linewidth',1.2);
             %     ylabel('Celcius')
             %     set(gca,'fontsize',14);
             xrange_p1 = get(gca,'xlim');
             yrange = get(gca,'ylim');
             ylim([yrange(1), yrange(1)+1.2]);
             grid on
             title(DN)
             %hold off
             
             %%====== Scalogram ======== %%
             subplot(2,2,3)
             levs = max(mag_wt(:)).*[0.5:0.1:1];
             pcolor(dist_equal, wvlen , mag_wt);
             shading flat
             hold on
             plot(DD(cond), WL(cond),'om','MarkerSize',10);
             [c,h] = contour(dist_equal, wvlen, mag_wt, levs,'k');
             axis tight
             xlabel('km');
             ylabel('wavelengths (km)');
             set(gca,'yscale','linear');
             hb = colorbar;
             set(hb,'location','southout','orientation','horizontal');
             set(get(hb,'xlabel'),'String','magnitude');
             set(gca,'ytick',[5, 10, 50,100],'TickDir','both')
             grid on
             xlim(xrange_p1 )
             set(gca,'fontsize',14);
             title('Magnitude Scalogram');
             hold off
             
         end
         
         % pack output data in a structure:
         prepped_data.SST_lowpassed = SST_lowpassed;
         prepped_data.moving_flag = moving_flag;
         prepped_data.dist_equal = dist_equal;                 % equally spaced distance axis
         prepped_data.pks = pks;                               % peaks in the SST record
         prepped_data.trghs = trghs;                           % troughs in the SST record
         prepped_data.mag_wt = mag_wt;                         % magnitude of wavelet transformation
         prepped_data.WL = WL;                            
         prepped_data.traj = traj;
         prepped_data.SST_lowpassed_dtr = SST_lowpassed_dtr;
         prepped_data.U10N_interp = U10_interp;
         prepped_data.U10N_lowpassed = U10_lowpassed;
         
         data.U10N_interp = U10_interp;
         data.U10N_lowpassed = U10_lowpassed;
         data.SST_interp = SST_interp;
         data.dist_equal = dist_equal;
         obj.Value = data;
         
         
         
            function data_lowpassed = lowpass_filter(data, wvlen_cutoff, spatial_res)
                % input data is expected to be equally spaced.
                %%% -- detrend:
                data_interp_dtr = detrend(data,1,'omitnan');                                % detrend (I don't need this.)
                trend_linear_ = data - data_interp_dtr;
                
                
                %%% -- low-pass filter:
                % wpass: normalized passband frequency: pi rad/sample.
                % filter out wavelength smaller than 20km;
                % wavenumber: 2*pi rad/20km;  resolution: 2km per sample
                % rad per sample = wavenumber_cutoff * resolution
                wn_cutoff_ = 2/wvlen_cutoff ;  % pi rad/km
                wpass_ = wn_cutoff_ * spatial_res;
                data_lowpassed_dtr = lowpass(data_interp_dtr, wpass_);   % retain frequency lower than the cutoff;
                data_lowpassed = data_lowpassed_dtr+trend_linear_;
                
                
            end
         
        end
        
        
        %%%% ========= function: select data segments for SST-Wind correlation analysis.  ========= %%%%
        % Purpose: 
        function [data_segs] = select_straight_traj_segments(obj, crits, sstname )
            
            %
            data = obj.Value;
            
            % 1. direction change <10°
            trajdir_cart= data.trajdir_cart;
            trajdir_cart(trajdir_cart<0) = trajdir_cart(trajdir_cart<0)+360;
            
            crit1 = true(size(trajdir_cart));
            crit_tmp = cosd(diff(trajdir_cart))>=cosd(crits.ddir);
            crit1(2:end) = crit_tmp;
            
            
            % first break the data down to individual segments:
            % do the segment separation with different along-wind,
            % cross-wind criteria:
            
            
            % locations of the data that satistied the criteria.
            locIDs = find(crit1);
            
            % find the edge of different non-continuous segments:
            edgeIDs = find(diff(locIDs)~=1);
            if ~iscolumn(edgeIDs)
                edgeIDs = edgeIDs';
            end
            
             fieldn = fieldnames(data);
             if ~isempty(locIDs)
                 if ~isempty(edgeIDs)
                     dimDist = [edgeIDs(1); diff(edgeIDs); length(locIDs) - edgeIDs(end)];
                     
                 else
                     dimDist = [length(locIDs)];
                 end
                 for i = 1:length(fieldn)
                     FN = fieldn{i};
                     val = data.(FN);
                     
                     if iscolumn(data.traj)
                         tmp = mat2cell(val(locIDs),dimDist, 1);
                     else
                         tmp = mat2cell(val(locIDs), 1, dimDist);
                     end
                     
                     for k = 1:length(tmp)
                         data_tmp(k).(FN) = tmp{k};
                     end
                     
                 end
             end
            
            figure
            subplot(3,1,1)
            
            plot(data.traj, trajdir_cart,'.b');
            hold on;
            plot(data.traj, trajdir_cart,'-k');
            set(gca,'fontsize',14);
            
            % for each segment, label it to be either a along-wind segment
            % or a cross wind segement based on the misalignment angle.
            % Note: a segment can not be both type.
            
            % 2. misalignment angle between the wind and the trajectory:
            % wind direction is in the atmospheric convention:
            wind_direction_name_list = {'wind_direction','wdir'};
            wdir_name_bool = ismember(wind_direction_name_list, fieldn);
            wdir_VN = wind_direction_name_list{wdir_name_bool};
            
            seg_label = zeros(length(dimDist),1);
            for k = 1:length(dimDist)
                trajdir_cart= data_tmp(k).trajdir_cart;
                trajdir_cart(trajdir_cart<0) = trajdir_cart(trajdir_cart<0)+360;
                
                winddir_cart = get_cartesian_direction(data_tmp(k).(wdir_VN),'Meteo');
                misang = abs(winddir_cart - trajdir_cart);
                
                plot(data_tmp(k).traj, misang,'-c');
                plot(data_tmp(k).traj, ones(size(data_tmp(k).traj)).*(crits.misang),'--k');
                plot(data_tmp(k).traj, ones(size(data_tmp(k).traj)).*(180-crits.misang),'--k');
                
                crit2 = (sind(misang)<=sind(crits.misang));   % 0~+/-crits.misang --> along_wind
                
                mask.along_wind =  crit2;
                mask.crx_wind = ~crit2;
                
                along_wind_rec = numel(find(mask.along_wind));
                crx_wind_rec = numel(find(mask.crx_wind));
                
                if along_wind_rec>crx_wind_rec
                    % label segment as alongwind
                    seg_label(k) = 1;
                    plot(data_tmp(k).traj, trajdir_cart,'+r');
                else
                    seg_label(k) = 0;
                    plot(data_tmp(k).traj, trajdir_cart,'*g');
                end
                
                set(gca,'ytick',[0:90:360]);
                grid on
                
                % compute segment length:
                data_tmp(k).segment_length = data_tmp(k).traj(end) - data_tmp(k).traj(1);
            end
            title([obj.ATOMIC_platform], [datestr(data_tmp(1).time(1)) ' ~ ' datestr(data_tmp(end).time(end))]);
            
            % separate segments into two cateogories based on the segment
            % labels;
            % in addition, drop the shorter segments.
            crit3 = [data_tmp.segment_length]>=crits.length;       % longer than 150 km
            crit3 = crit3';
            mask.along_wind = (seg_label==1)&crit3;
            mask.crx_wind = (seg_label==0);                       % the distance criterion on applies to the along-wind segments.
            
            
            seg_type = fieldnames(mask);
            for s = 1:length(seg_type)
                ST = seg_type{s};
                data_segs.(ST) = data_tmp(mask.(ST));
                
                % sanity check
                subplot(3,1,s+1);
                if ~isempty(data_segs.(ST))
                for k = 1:length(data_segs.(ST))
                    hold on
                    plot(data.lon, data.lat,'-k');

                    scatter(data_segs.(ST)(k).lon, data_segs.(ST)(k).lat, 20, ...
                        data_segs.(ST)(k).(sstname),'filled');
                    theta = get_cartesian_direction(data_segs.(ST)(k).(wdir_VN),'Meteo');
                    u = 0.01 .* cosd(theta);
                    v = 0.01.* sind(theta);
                    quiver(data_segs.(ST)(k).lon(1:4:end), data_segs.(ST)(k).lat(1:4:end),...
                        u(1:4:end),v(1:4:end),'k');
                    text(mean(data_segs.(ST)(k).lon), mean(data_segs.(ST)(k).lat)-0.02, ...
                        ['seg' num2str(k)],'fontsize',12);
                    plot(data_segs.(ST)(k).lon(1), data_segs.(ST)(k).lat(1),'*m','markersize',13);
                end
                axis('equal');
                if s == 1
                title(['along wind segments = ' num2str(length(data_segs.(ST)))], ...
                    [datestr(data_segs.(ST)(1).time(1)) ' ~ ' datestr(data_segs.(ST)(end).time(end))]);
                end
                set(gca,'fontsize',14);
                
                end
                
            end
        end
        
        %%%%% =========== function: re-orient data sequence so it points to 
        %%%%% ===========            downwind direction.       ======= %%%%
        function obj_new = reorient_data_sequence(obj)
            % determine wheter or not the data sequence needs to be
            % rotated:
            data = obj.Value;
            obj_new = obj;
            
            data_new = data;
            
            nsegs = length(data);
            
            wind_direction_name_list = {'wind_direction','wdir'};
            fieldn = fieldnames(data);
            wdir_name_bool = ismember(wind_direction_name_list, fieldn);
            wdir_VN = wind_direction_name_list{wdir_name_bool};

            
            for iseg = 1:nsegs
                
                trajdir_cart = data(iseg).trajdir_cart;
                winddir_cart = get_cartesian_direction(data(iseg).(wdir_VN),'Meteo');
                
                % compare the trajectory direction with the wind direction;
                if cosd(trajdir_cart-winddir_cart)<0
                    % rotate all the variables in the data object, except:
                    % lat, lon, traj, time, local_time, trajdir_cart(what else?)
                    excp_fields = {'lon','lat','traj','time','local_time','distance_axis'};
                    fieldn = fieldnames(data(iseg));
                    
                    [ny, nx] = size(data(iseg).lon);
                    
                    for i = 1:length(fieldn)
                        x = fieldn{i};
                        if ~ismember(x, excp_fields)
                            dimcrit = size(data(iseg).(x)) == size(data(iseg).lon);
                            if dimcrit
                                if iscolumn(data(iseg).(x))
                                    data_new(iseg).(x) = flipud(data(iseg).(x));
                                else
                                    data_new(iseg).(x) = fliplr(data(iseg).(x));
                                end
                                
                                figure(12);
                                if strcmp(x, 'sea_water_temperature') || strcmp(x,'tsea') || strcmp(x,'wspd_10N')
                                    plot(data(iseg).(x));

                                    hold on;
                                    plot(data_new(iseg).(x),'r');

                                    hold off;
                                pause
                                end
                                
                            else
                                % find out which dimension is the same;
                                dimID = find(dimcrit==1);
                                if dimID ==1
                                    data_new(iseg).(x) = flipud(data(iseg).(x));
                                else
                                    data_new(iseg).(x) = fliplr(data(iseg).(x));
                                end
                                
                            end
                        end
                    end
                    disp(['seg#', num2str(iseg,'%2.2i') 'sequence re-oriented.']);

                    
                else
                    disp('no need for re-orientation.');
                end
            end
            
            obj_new.Value = data_new;
        end
        
        %%%% =================   function #2b   ===================== %%%%
        %%%% Purpose: project data onto the RHB or wave glider trajectory.
        function  obj = project_wind_onto_trajectory(obj, opt)
            switch nargin
                
                case 1
                    copt = false;
                    %wopt = true;
                case 2
                    copt = opt;
                   % wopt = true;
                
            end
            
            
            data = obj.Value;
            
            nsegs = length(data);
            
            
            wind_direction_name_list = {'wind_direction','wdir'};
            wind_speed_namelist = {'wind_speed','wspd_10N'};
            
            fieldn = fieldnames(data(1));
            wdir_name_bool = ismember(wind_direction_name_list, fieldn);
            wdir_VN = wind_direction_name_list{wdir_name_bool};
            
            wspd_name_bool = ismember(wind_speed_namelist, fieldn);
            wspd_VN = wind_speed_namelist{wspd_name_bool};
            
            for iseg = 1:nsegs
                [wu_alg, wv_crx] = compute_algtraj_crxtraj_velocity_component(wspd_VN, wdir_VN);
                data(iseg).u_algtraj = wu_alg;
                data(iseg).v_crxtraj = wv_crx;
            end
            
            
            if copt   % project current as well:  % no current from the wave glider.
                current_speed_namelist = {'cspd','drift_speed'};
                current_direction_namelist = {'cdir','drift_direction'};
                
                fieldn = fieldnames(data(1));
                cdir_name_bool = ismember(current_direction_namelist, fieldn);
                cdir_VN = current_direction_namelist{cdir_name_bool};
                
                cspd_name_bool = ismember(current_speed_namelist, fieldn);
                cspd_VN = current_speed_namelist{cspd_name_bool};
                
                for iseg = 1:nsegs
                    [cu_alg, cv_crx]=compute_algtraj_crxtraj_velocity_component(cspd_VN, cdir_VN);
                    data(iseg).cu_algtraj = cu_alg;
                    data(iseg).cv_crxtraj = cv_crx;
                end
            end
                
            obj.Value = data;
            
           %% turn the following code into a nested function:           
           % for iseg = 1:nsegs
            function [u_algtraj, v_crxtraj]=compute_algtraj_crxtraj_velocity_component(spd_VN, dir_VN)   
                trajdir_cart = data(iseg).trajdir_cart;
                trajdir_cart(trajdir_cart<0) = trajdir_cart(trajdir_cart<0)+360;

                % although the variable name is called "wind direction" it
                % also stores the current direction information in it. 
                dir_cart = get_cartesian_direction(data(iseg).(dir_VN),'Meteo');
                
                dtheta = trajdir_cart - dir_cart;
                against_locs = cosd(dtheta)<0;
             
                ang = dtheta;
                %ang(against_locs) = (180-abs(dtheta(against_locs))).*sign(dtheta(against_locs));  % I doubt this again.
                ang(against_locs) = (dtheta(against_locs));                % updated to this version (Jan 24, 2022, the wind direction of the project along traj seems correct this way.
                
                u_algtraj = data(iseg).(spd_VN) .* cosd(ang);              % <0
                v_crxtraj = data(iseg).(spd_VN) .* sind(ang);              % > 0?
                
                % sanity check;
                utmp = data(iseg).(spd_VN) .* cosd(dir_cart);
                vtmp = data(iseg).(spd_VN).* sind(dir_cart);
                
                ualg_x = u_algtraj .* (cosd(trajdir_cart));                % is the direction of this correct??
                ualg_y = u_algtraj .* sind(trajdir_cart);
                
                ucrx_x = v_crxtraj .* (cosd(trajdir_cart-90));
                ucrx_y = v_crxtraj .* sind(trajdir_cart-90);
                
               

                
                figure(10);clf;
                quiver(data(iseg).lon(1:5:end), data(iseg).lat(1:5:end), ...
                    utmp(1:5:end), vtmp(1:5:end),'k');
                hold on;
                quiver(data(iseg).lon(1:5:end), data(iseg).lat(1:5:end), ...
                    ualg_x(1:5:end), ualg_y(1:5:end),'r');
                hold on;
                quiver(data(iseg).lon(1:5:end), data(iseg).lat(1:5:end), ...
                    ucrx_x(1:5:end), ucrx_y(1:5:end),'c');
                title({[spd_VN '; iseg=', num2str(iseg)];[datestr(data(iseg).time(1))]});
                axis('equal');
                pause
            end
            
            
            
        end
        
        %%%% =================   function # ??   ===================== %%%%
        %%%% Purpose: compute the along-traj Rossby number from the surface
        %%%% current measurements. Ro = eta/f; eta_algtraj = dv_crstraj/dx 
        %%%% x: along trajectory direction, y: cross-trajectory direction.
        %%%% Date: Jan 24, 2022
        
        function obj_new = compute_along_traj_current_relative_vorticity_and_Ro(obj)
            
            % note: input obj only has 1 segment. 
            % project velocity on to the trajectory and find the
            % corresponding u_algtraj, v_crxtraj component
            % better to call a function for this purpose.
            omg = 2*pi/86400;   % rad/s
            curflag = true;
            obj = project_wind_onto_trajectory(obj, curflag);
            
            obj_new = obj;

            data = obj.Value;
            nsegs = length(data);
            
            for iseg = 1:nsegs
                distance = data(iseg).distance_axis;
                
                % take out the current velocity:
                %cur_alg = data.cu_algtraj;
                cur_crx = data(iseg).cv_crxtraj;
                
                
                % compute the vertical vorticity along trajectory
                % finite difference, (the vorticitiiy is computed at the
                % mid-point of the distance axis.)
                dist_mid = 0.5*(distance(1:end-1)+distance(2:end));
                eta_algtraj_mid = diff(cur_crx)./(diff(distance)*1E3);
                
                % interpolate back to the same grid as the "traj" variable.
                eta_algtraj = interp1(dist_mid, eta_algtraj_mid, distance);
                
                
                % normalize the vertical vorticity by the Coriolis
                % parameter:
                f = 2*omg*sind(data(iseg).lat);
                Ro = eta_algtraj./f;
                data(iseg).vort_algtraj = eta_algtraj;
                data(iseg).Ro = Ro;
            end
            
           obj_new.Value = data;

                
            
        end
        
        
        %%%% =================   function # ??   ===================== %%%%
        %%%% Purpose: break long data down into sections
%         function obj_p = break_record_into_pieces(obj, l)
%             % input: l: length of each pieces
%             for i = 1:length(obj)
%                 data = obj(i).Value;
%                 
%                 seglen = data.segment_length;
%                 
%                 ndiv = 2;   % starting with this basic dividing number
%                 operation_flag = true;
%                 
%                 if seglen>150
%                     % require operation
%                     while operation_flag
%                         
%                     end
%                     
%                 else
%                    % keep as is and count as a segment.
%                    obj_p(cnt).data = ;
%                 end
%                 
%             end
%             
%         end
       

        %%%%% =========== Utility function: saving an output;  ======= %%%%
        function save_data(obj, svdataInfo)
           absFN = [svdataInfo.svdir filesep svdataInfo.filename];      
            feval('save', absFN, svdataInfo.varname);
        end


        
    end
    
    
end