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
    end

    % establish method modules:
    methods
        
        % establish functions:
        %%%% ===================  function #0  ===================== %%%%
        %%%% initialize the object:
        function obj = ATOMIC_Wavelet(ts1, ts2, traj,time, ATOMIC_platform)
            if nargin > 0 
               obj.ts1 = ts1;
               obj.ts2 = ts2;
               obj.traj = traj;
               obj.time = time;
               obj.ATOMIC_platform = ATOMIC_platform;
         
            else
                disp('initialize object with idealized sin/cos waves');
                theta = -5*pi:0.01*pi:5*pi;
                noise = randn(size(theta));
                
                obj.ts1 = sin(theta)+noise;
                obj.ts2 = cos(theta)+noise;
                obj.traj = 1:length(theta);
                obj.time = linspace(0,60, length(theta));
                obj.ATOMIC_platform = 'none';
            end
        end
        
        %%%% ===================  function #1  ===================== %%%%
        %%%% Purpose: plot the original input records;
        function hfigs= plot_data(obj,labelstr)
            x = obj.ts1;
            y = obj.ts2;
            traj = obj.traj;
            time = obj.time;
            
            x_dtr = detrend(x,1);
            y_dtr = detrend(y,1);

            
            % purpose: plot the two input "time series" for a view.
            hfigs = figure;
            subplot(2,1,1)
            yyaxis left
            plot(traj,x,'linewidth',1.2);
            ylabel(labelstr.ts1);
            
            yyaxis right
            plot(traj,y,'linewidth',1.2);
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
        function [hfigs, wtcout]= wavelet_coherence_toolbox(obj)
            % output key parameters from wavelet coherence.  
            %[Rsq,period,scale,coi,sig95]=wtc(x,y,[,settings])
            addpath '/Users/xchen/Documents/MATLAB/customized_functions/wavelet-coherence-master';
            
            x = obj.ts1;
            y = obj.ts2;
            traj = obj.traj;
            time = obj.time;
            
            % the input x, and y has been flipped to progress downwind;
            x_dtr = detrend(x,1);
            y_dtr = detrend(y,1); 
            
           [Rsq, period, scale, coi, wtcsig]=wtc(x_dtr, y_dtr);           % figure out what do period and scale stand for. 

            hfigs = figure;
            wtc(x_dtr, y_dtr);
            ax = gca;
            xtick = ax.XTick;
            ax.XTickLabel = xtick.*2;
            ytick = ax.YTick;
            ax.YTickLabel = 2*2.^ytick;
            ax.YDir = 'normal';
            xlabel('Downwind Distance (km)');
            ylabel('Wavelength (km)');
            set(gca,'fontsize',14)
            ax.Title.String=[obj.ATOMIC_platform];
            ax.TickDir = 'both';
            
            
            wtcout.coherence_squared = Rsq; 
            wtcout.period = period;       % think about a better name..
            wtcout.scale = scale;
            wtcout.coi = coi;
            wtcout.sig95 = wtcsig; 
            
            % need to output phase angle as well:
            % wtcout.phase = phase_ang;
            
%             varargout={Rsq,period,scale,coi,wtcsig};
%             varargout=varargout(1:nargout);
        end
        
        
        %%%% ===================  function #4  ===================== %%%%
        %%%% Purpose: obtain the scale and phase in area with high and significant area;        
        function data_o=get_averaged_scale_and_phase_from_high_coherence_region(obj)
            %    - find local maximum in the coherence
            %    - conditional selection: select local coherence maximum that is higher
            %    than the sig95 threshold.
            %    - get the averaged phase angle and wavelength from the region with
            %    high coherence
            
            
        end
        
        
        
        
    end
    
end
