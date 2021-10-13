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
            tmin = floor(min(ds.time));
            tmax = ceil(max(ds.time));
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
        

            
        %%%% =================   function #2   ===================== %%%%
        %%%% Purpose: compute distance traveled by the ship/instrument
        function obj = compute_distance_travelled(obj)
            
            data = obj.Value;     % this is a matlab structure;
            
            dlat = diff(data.lat);
            dlon = diff(data.lon);
            
            % average latitude between to consecutive ship locations:
            mid_lat = 0.5*(data.lat(1:end-1) + data.lat(2:end));
            xdist = abs(dlon) .* 111.*cosd(mid_lat);  % units: km
            ydist = abs(dlat) .* 111;
            
            % construct trajectory coordinate, (0,0) at the first location of the
            % ship.
            traj_x = [0, cumsum(xdist)];
            traj_y = [0, cumsum(ydist)];
            traj = sqrt(traj_x.^2 + traj_y.^2);
            
            data.traj = traj; 
            
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
                fignum = hfig_check.Number;
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

                figname = [DN '_SSTgrad_segment_v0_lowpassed.jpg']
                xc_savefig(gcf, figsvdir, figname, [0 0 14 10]);
            end
            
        end
        
       
        %%% =================== function #6 ====================== %%%
        function obj_new = map_to_distance_axis(obj, spatial_res)
            data = obj.Value;
            
            spd_thres = 1;      % 1m/s-> 3.6km/hr;
            moving_flag = get_ship_speed_mask(obj, spd_thres);
            traj = data.traj(moving_flag);
            dist_equal = [0:spatial_res:round(max(traj))];
            
            % interpolation for variaous values:
            fieldn = fieldnames(data);
            for i = 1:length(fieldn)
                FN = fieldn{i};
                tmp_val = data.(FN)(moving_flag);
                data_interp.(FN) = interp1(traj, tmp_val, dist_equal);
            end
            data_interp.distance_axis = dist_equal;
            
            
            
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
         U10_interp = interp1(traj, data.wspd_10N(moving_flag),  dist_equal, 'linear');
         
         % the wind direction is likely not sensitive to the speed of the
         % ship(?)
         data.ship_wind_angle_interp  = interp1(data.xmid, data.ship_wind_angle, dist_equal);
         
         
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
                data_interp_dtr = detrend(data,1);                                % detrend (I don't need this.)
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
        
        
        
        
        %%%%% =========== Utility function: saving an output;  ======= %%%%
        function save_data(obj, svdataInfo)
           absFN = [svdataInfo.svdir filesep svdataInfo.filename];      
            feval('save', absFN, svdataInfo.varname);
        end


        
    end
    
    
end