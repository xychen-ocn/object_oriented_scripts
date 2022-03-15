% ======================================================================= %
% class: ATOMIC_Wavelet
% purpose: Wavelet workflow (produce results related to wavelet coherence)
% Author: Xuanyu Chen 
% Date: v0 @ 10/15/2021
% ======================================================================= %

classdef ATOMIC_Wavelet
    properties
        ts1 
        ts2
        traj
        time
        FigOutput_rootdir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/object_oriented_scripts/Figs';
        ATOMIC_platform
        xres
        wtcstat
    end

    % establish method modules:
    methods
        
        % establish functions:
        %%%% ===================  function #0  ===================== %%%%
        %%%% initialize the object:
        function obj = ATOMIC_Wavelet(ts1, ts2, traj,time, ATOMIC_platform, xres)
            if nargin > 0 
               obj.ts1 = ts1;
               obj.ts2 = ts2;
               obj.traj = traj;
               obj.time = time;
               obj.ATOMIC_platform = ATOMIC_platform;
               obj.xres = xres;
         
            else
                disp('initialize object with idealized sin/cos waves');
                %theta = -5*pi:0.01*pi:5*pi;
                x = -10:0.1:10;
                wn = 2*pi/20;
                theta = wn.*x;
                %noise = randn(size(theta));
                noise = 0;
                obj.ts1 = sin(theta)+noise;
                obj.ts2 = cos(theta)+noise;
                obj.traj = 1:length(theta);
                obj.time = datenum(now)+[0:length(theta)]/24;
                obj.ATOMIC_platform = 'none';
                obj.xres = 0.1;
            end
        end
        
        %%%% ===================  function #1  ===================== %%%%
        %%%% Purpose: plot the original input records;
        function hfigs= plot_data(obj,labelstr)
            x = obj.ts1;
            y = obj.ts2;
            traj = obj.traj;
            time = obj.time;
            
%             x_dtr = detrend(x,0);
%             y_dtr = detrend(y,0);
            
            if strcmp(obj.ATOMIC_platform, 'none')
                x_dtr = detrend(x,0);
                y_dtr = detrend(y,0);
            else
                
                x_dtr = detrend(x,1);
                y_dtr = detrend(y,1);
            end

            
            % purpose: plot the two input "time series" for a view.
            hfigs = figure(1); clf;
            subplot(2,1,1)
            yyaxis left
            h(1) = plot(traj,x,'linewidth',1.2);
            hold on
            plot(traj, x-x_dtr,'--','color',h(1).Color);
            ylabel(labelstr.ts1);
            
            yyaxis right
            h(2) = plot(traj,y,'linewidth',1.2);
            hold on
            plot(traj, y-y_dtr,'--','color',h(2).Color);
            xlabel('downwind distance (km)');
            title([datestr(min(time)) '~' datestr(max(time)) ('UTC')])
            ylabel(labelstr.ts2);
            set(gca,'fontsize',14);
            grid on
            
            
            subplot(2,1,2)
            yyaxis left
            plot(traj,x_dtr,'linewidth',1.2);
            ylabel(['\Delta' labelstr.ts1]);
            %ylim([-0.5, 0.5]);
            grid on
            
            yyaxis right
            plot(traj,y_dtr,'linewidth',1.2);
            xlabel('downwind distance (km)');
            ylabel(['\Delta' labelstr.ts2{end}]);
            
            title('Anomalies');
            set(gca,'fontsize',14);
            
            
        end
 
       
        %%%% ===================  function #2  ===================== %%%%
        %%%% Purpose: transform input records to normal distribution. (via
        %%%% box-cox power transformation)
        function [obj_trans, hfig] = transform_to_normal_distribution(obj, labelstr)
            [trans.ts1, lmb_x, maxlog_x] = boxcox(obj.ts1);
            [trans.ts2, lmb_y, maxlog_y] = boxcox(obj.ts2);
            
            obj_trans = obj;
            obj_trans.ts1 = trans.ts1;
            obj_trans.ts2 = trans.ts2;
            
            traj = obj.traj;
            
            % ideailly, I should test if it is a normal distribution.
            
            % add plots to compare the time series before transformation
            % and after transformation;
            varn = {'ts1','ts2'};
            
            hfig = figure(10);clf;
            for i = 1:2
               % titlestr = strsplit(labelstr.(varn{i}),'');
                
                subplot(2,1,i);
                yyaxis left
                plot(traj, obj.(varn{i}),'linewidth',1.2);
                ylabel(labelstr.(varn{i}));
                %ylim([25.8, 27.2]);
                
                yyaxis right               
                plot(traj, trans.(varn{i}),'linewidth',1.2);
                ylabel('transformed data');
                %title( titlestr{i} ) 
                xlabel('Downwind distance (km)');                
                set(gca,'fontsize',14);
               
                
                axis_parent = gca;
                pos_parent = get(gca,'pos');
                pos_ax1 = [pos_parent(1)+0.6* pos_parent(3), pos_parent(2)+0.05*pos_parent(4), pos_parent(3)*0.15, pos_parent(4)*0.25];
                ax1 = axes('pos',pos_ax1);
                
                % plot distribution on sub axes in the same panel.
                histogram(ax1, obj.(varn{i}),6);
                title('original data');
                
                pos_ax2 = [pos_parent(1)+0.8*pos_parent(3), pos_parent(2)+0.05*pos_parent(4), pos_parent(3)*0.15, pos_parent(4)*0.25];
                ax2 = axes('pos', pos_ax2);
                histogram(ax2, trans.(varn{i}),6)
                title('box-cox transformed');
                set(gca,'fontsize',12);
            end
            
            
            
        end
        
        
        %% TO_UPDATE!!! (The function below is for wavelet coherence and need to be updated to produce a 
        % bunch of figures for our interest.
        %%%% ===================  function #3  ===================== %%%%
        %%%% Purpose: use the wavelet toolbox by Grinsted et al. (2004) to
        %%%% compute wavelet coherence. 
        function [hfigs, wtcout, outCOI_stat, inCOI_stat]= wavelet_coherence_toolbox(obj)
            % output key parameters from wavelet coherence.  
            %[Rsq,period,scale,coi,sig95]=wtc(x,y,[,settings])
            addpath '/Users/xchen/Documents/MATLAB/customized_functions/wavelet-coherence-master';
            
            x = obj.ts1;
            y = obj.ts2;
            traj = obj.traj;
            time = obj.time;
            
            % the input x, and y has been flipped to progress downwind;
            if strcmp(obj.ATOMIC_platform, 'none')
                x_dtr = detrend(x,0);
                y_dtr = detrend(y,0);
            else
                
                x_dtr = detrend(x,1);
                y_dtr = detrend(y,1);
            end
%             if iscolumn(x)
%                 ds1=[traj, x_dtr];
%                 ds2=[traj, y_dtr];
%             else
%                 ds1=[traj', x_dtr'];
%                 ds2=[traj', y_dtr'];
%             end
            
           [Rsq, normalized_period, normalized_scale, coi, wtcsig, Wxy]=wtc(x_dtr, y_dtr);           % figure out what do period and scale stand for. 

            hfigs = figure(11); clf;
            % set size of the figure;
            set(gcf, 'pos',[616   495   653   523]);
            wtc(x_dtr, y_dtr);
            ax = gca;
            xlabel('Downwind Distance (km)');
            ylabel('Wavelength (km)');
            set(gca,'fontsize',14)
            
            pause(0.05);
            xtick = ax.XTick;
            
            ax.XTickLabel = xtick.* obj.xres;
            ytick = ax.YTick;
            ax.YTickLabel = 2.^ytick .* obj.xres;
            ax.YDir = 'normal';
                       
            ax.Title.String=[obj.ATOMIC_platform];
            ax.TickDir = 'both';
            
            
            wtcout.coherence_squared = Rsq; 
            wtcout.period = normalized_period;       % units: sample
            wtcout.wavelet_scale = normalized_scale;          % S0*2^(j*DJ), S0: the smallest scale of the wavelet  = 2*DT or 2*DX (in this case), DJ: spacing between discrete scales, default is 0.25
            wtcout.coi = coi;
            wtcout.sig95 = wtcsig; 
            wtcout.Wxy = Wxy;                                 % wavelet cross-spectrum
            
            % call a function to get the averaged scale and phase from the
            % high coherence region.
            [outCOI_stat, inCOI_stat] = obj.get_averaged_scale_and_phase_from_high_coherence_region(wtcout, 0.4);
            ds_out.wvlet_coh = wtcout;
            ds_out.outCOI_metrics = outCOI_stat;
            ds_out.inCOI_metrics = inCOI_stat;
            obj.wtcstat = ds_out;
        end
        
        
        %%%% ===================  function #4  ===================== %%%%
        %%%% Purpose: obtain the scale and phase in area with high and significant area;        
        function [notCOI, COI_area]=get_averaged_scale_and_phase_from_high_coherence_region(obj, wtc_stat, varargin)
            %    - find local maximum in the coherence
            %    - conditional selection: select local coherence maximum that is higher
            %    than the sig95 threshold.
            %    - get the averaged phase angle and wavelength from the region with
            %    high coherence
            
            if nargin > 1
                thres = varargin{1};
            else
                
                thres = 0.6;     % this threshould number control whether or not to count the signal influenced by the edge as
                % a valid signal. Can be varied to relax
                % criteria.
            end
            disp(['thres=' num2str(thres)]);
            x = obj.ts1;
            y = obj.ts2;
            traj = obj.traj;
            time = obj.time;
            xres = obj.xres;
            
            [xx, yy] = meshgrid(traj, log2(wtc_stat.period));
            [tt, ~] = meshgrid(time, log2(wtc_stat.period));
            
            Rsq = wtc_stat.coherence_squared;
            sig95_locs = wtc_stat.sig95 >=1;  % the definition of sig95 was wtcsig = Rsq./wtcsig in the code..   
            
            figure(11);
            pcolor(xx,yy, double(sig95_locs)); shading flat; hold on;
            
            % obtain averaged phase and scale within the regions of interests:
            % xv, yv represent the contour (polygon) of the 95% significant
            % region. xx, yy is the two axis I had.
            
            % use regionprop to find the area with 95% significance:
            CCr = bwconncomp(sig95_locs,8);
            Lr = labelmatrix(CCr);
            stats=regionprops(Lr, sig95_locs,'PixelIdxList','Image','Area');
            num_sig95_patches = length(stats);
            
            % use pixel index to compute the averaged phase and ...
            xcoi = [traj([1 1])-xres*.5, traj, traj([end end])+xres*.5];
            ycoi = log2([wtc_stat.period([end 1]) ...
                wtc_stat.coi ...
                wtc_stat.period([1 end])]);
            
            plot(xcoi, ycoi,'--w'); hold off;
            
            wvlen_scale = wtc_stat.period.*xres;
            [~, WVLEN] = meshgrid(traj, wvlen_scale);
            Wxy  = wtc_stat.Wxy; 
            
            cnt1=0;
            cnt2=0;
            for ii = 1:num_sig95_patches
                idx = stats(ii).PixelIdxList;
                xsel = xx(idx);
                ysel = yy(idx);
                tsel = tt(idx);                   % selected time window
                
                [inCOI, onCOI] = inpolygon(xsel, ysel, xcoi, ycoi);
                
                cond = ~inCOI & ~onCOI ;
                
                num_pixel_outCOI = length(idx(cond));
                
                good_frac = num_pixel_outCOI/length(cond);
                disp(num2str(good_frac));
                if round(good_frac*10)/10>=thres
                    % a good patch
                    % significant 95% location is outside of the COI;
                    cnt1 =cnt1+1;
                    notCOI.ave_Wxy(cnt1) = mean(Wxy(idx(cond)),'omitnan');
                    notCOI.ave_phase(cnt1) = angle(notCOI.ave_Wxy(cnt1));
                    notCOI.ave_wvlen(cnt1) = mean(WVLEN(idx(cond)));
                    notCOI.ave_cohsq(cnt1) = mean(Rsq(idx(cond)));
                    notCOI.area(cnt1) = stats(ii).Area * good_frac;
                    
                    % figure out the duration of the correlated signal.
                    [xloc.x0, stid] = min(xsel(cond));
                    [xloc.xN, edid] = max(xsel(cond));
                    
                    notCOI.duration(cnt1) = (tsel(edid)-tsel(stid))*24;    % units: hour
                    notCOI.physical_loc(cnt1) = xloc;
                    
                    figure(21);
                    polarplot(notCOI.ave_phase(cnt1), notCOI.ave_wvlen(cnt1) , ...
                             '.b','MarkerSize',notCOI.area(cnt1)*0.1);
                    hold on;
                    
                else
                    % a bad patch (compute the metrics from all pixels, even if those are outside COI.)
                    cnt2=cnt2+1;
                    COI_area.ave_Wxy(cnt2) = mean(Wxy(idx), 'omitnan');
                    COI_area.ave_phase(cnt2) = angle( COI_area.ave_Wxy(cnt2));
                    COI_area.ave_wvlen(cnt2) = mean(WVLEN(idx));
                    COI_area.ave_cohsq(cnt2) = mean(Rsq(idx));
                    COI_area.fraction(cnt2) = good_frac;                   
                    COI_area.area(cnt2) = stats(ii).Area * (1-good_frac);
                   
                   
                end
            end
            
            if cnt1==0
                notCOI.ave_phase=NaN;
                notCOI.ave_wvlen=NaN;
                notCOI.ave_cohsq=NaN;
                notCOI.ave_Wxy=NaN;
                notCOI.area = NaN;
                notCOI.duration = NaN;
                notCOI.physical_loc.x0 = NaN;
                notCOI.physical_loc.xN = NaN;
            end
            
            if cnt2==0
                COI_area.ave_phase=NaN;
                COI_area.ave_wvlen=NaN;
                COI_area.ave_cohsq=NaN;
                COI_area.fraction =NaN;
                COI_area.ave_Wxy = NaN;
                COI_area.area =NaN;
            end
            
            
            %
        end
        
        %%%% ===================  function #4  ===================== %%%%
        %%%% Purpose: make scatter plots of scale and phase in area with
        %%%% high and significant coherence in a polar coordinate.        
        function [hfig, P90] = plot_scale_and_phase(obj,pthres, varargin)
            
           % enalbe optional input argument:
           num_optargs = length(varargin);
          
           switch num_optargs
               case 0
                   inCOI_flag = false;
                   colorcode_flag = true;
               case 1
                   inCOI_flag = varargin{1};
                   colorcode_flag = true;
               case 2
                   inCOI_flag = varagin{1};
                   colorcode_flag = varargin{2};
                   
                  
           end
           
           %  get the data metrics out of object:
           outCOI_stat = obj.wtcstat.outCOI_metrics;  %not influenced by the edge effects;
           inCOI_stat = obj.wtcstat.inCOI_metrics;
           
           
           wvlens.QC = [outCOI_stat.ave_wvlen];
           phase.QC = [outCOI_stat.ave_phase];
           [areasize.QC, sid] = sort([outCOI_stat.area],'descend');
           duration.QC = [outCOI_stat.duration];
           totalnum_locs.QC = length(wvlens.QC);
           
           wvlens.QC=wvlens.QC(sid);
           phase.QC=phase.QC(sid);
           duration.QC=duration.QC(sid);
           
           if inCOI_flag
               wvlens.bad = [inCOI_stat.ave_wvlen];
               phase.bad = [inCOI_stat.ave_phase];
               areasize.bad = [inCOI_stat.area];
               %duration.bad = [inCOI_stat.duration];
               totalnum_locs.bad = length(wvlens.bad);
           end
           
           
           % compute the 90th percentile of wavelenth:
           disp(num2str(pthres));
           p90 = prctile(wvlens.QC, pthres);
           if inCOI_flag
               p90_all = prctile([wvlens.QC, wvlens.bad], pthres);
           end
           P90.QC = p90;
           P90.all = p90_all;
           disp(['p90=' num2str(p90)]);
           ocrit= wvlens.QC>p90;
           outliers.WL = wvlens.QC(ocrit);
           outliers.theta = phase.QC(ocrit);
          
            hsc=[];
            hfig =figure; 
             
            if colorcode_flag
                hsc(1)= polarscatter(phase.QC, wvlens.QC, sqrt(areasize.QC)*5, duration.QC, 'filled');
            else
                hsc(1)= polarscatter(phase.QC, wvlens.QC, sqrt(areasize.QC)*5, 'b', 'filled');               
            end
            hold on;
            
            if inCOI_flag
                hsc(2) = polarscatter(phase.bad, wvlens.bad, sqrt(areasize.bad)*5,'k', 'marker','x');
            end
            
           % plot the p90 radius:
            th = [0:0.1*pi:2*pi];            
            hl = polarplot(th, p90.*ones(size(th)),'--r','linewidth',2); 

            
            %ho = polarplot(outliers.theta, outliers.WL, 'ow','markersize', 10);
            colormap(turbo)
            hb = colorbar;
            set(get(hb,'xlabel'),'String','data duration (hours)');
            %hb.Title.String = '10^{\^}';
            %title(['wavelet coherence between SST and ' varn_long]);
            set(gca,'fontsize',14);
            rlim([0 40])
            hold off;
            if length(hsc)==1
                legend([hsc,  hl],{'Good data',['P' num2str(pthres,'%i') '^{th}=' num2str(round(p90),'%3.1f') 'km']},'fontsize',14);
            else
                legend([hsc,  hl],{'Good data','Edge Affected',['P' num2str(pthres,'%i') ...
                   '^{th}=' num2str(round(p90),'%3.1f') 'km']},'fontsize',14);

            end

            
    
            
        end
        
        
        
    end
    
end
