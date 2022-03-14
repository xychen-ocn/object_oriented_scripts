% recover scripts that I used to generate the noise or wavelet.
xres = 0.2;
x = 0:xres:200;
wn = 2*pi/20;
theta = wn.*x;
noise1 = randn(size(theta));
noise1_n = noise1./max(noise1);
noise2 = randn(size(theta));
noise2_n = noise2./(max(noise2));
phase =[0:pi/2:2*pi-pi/2];

G_env = ;   % Guassian envelope;

A = [3, 1];
x0 = 5+5*20;

figsvdir = 'Figs/demo';

k1 = wn;
k2 = 2*pi/100;
ts3 = A(2)*sin(k2.*x);
ts4 = A(2)*sin(theta)+ noise1_n*0.5+ts3;
plot(x,ts4)



for i = 1:length(phase)
    ph = phase(i);
    ts2 = A(2)*sin(theta)+ noise2_n*0.1; %+noise2_n;     % SST:
    ts1 = A(1)*sin(theta+ph)+noise1_n*0.5;     %wind
    ts2(x<(x0-0.75*20))=0;
    ts2(x>(x0+0.75*20))=0;
    ts2= ts2 + noise2_n;
    % ts2(x<35+20*phase/(2*pi))=-1;
    % ts2(x>55+20*phase/(2*pi))=-1;
    % ts1 = ts1+noise;
    % ts2 = ts2+noise;
    traj = x;
    time = datenum(now)+[0:length(theta)]/24;
    ATOMIC_platform = ['L=20, phase=' num2str(ph/pi) 'pi'];
    testobj = ATOMIC_Wavelet(ts1, ts2, traj, time, ATOMIC_platform, xres);
    labelstr.ts1 = {'sin(kx-\theta)', 'w/ random noise'};
    labelstr.ts2 = {'sin(kx)','w/ random noise'};
    testobj.plot_data(labelstr);
    xc_savefig(gcf, figsvdir,['signal_L20_phase' num2str(ph/pi) 'pi_new.jpg'],[0 0 8 6]);

    
    [hfig, wtc_stat(i), outCOI_stat(i), inCOI_stat(i)]=testobj.wavelet_coherence_toolbox;
    xc_savefig(hfig, figsvdir,['sine wave with random noise_L20_phase' num2str(ph/pi) 'pi_new.jpg'],[0 0 8 6]);
end


% think about creating a "monte carlo" simulation with these two simple
% signals, but changing the noise and location of the correlation, but keep
% the phase as the same. generate 100 case and then plot the results in
% polar plots.

%% create a polar plot to visualize the results:
figure();
th = [outCOI_stat.ave_phase];
s = [outCOI_stat.ave_wvlen];
msize =  [outCOI_stat.area];
for i = 1:length(msize)
    hp(i) = polarplot(th(i), s(i),'.b','MarkerSize',msize(i)*0.2);
    hold on;
end
set(gca,'fontsize',14);







%% generate 100 cases with different random noises and a randomly generated 
% phase angle and wavelength for the two signals.
%% now generate a time series with three bumps located at specified locations;
figsvdir = './create_demo/figs';
if ~exist(figsvdir)
    mkdir(figsvdir)
end
addpath('./create_demo');

x = 0:0.2:150;  % km
ncount = 100;
delta_phase = 10/180*pi;
delta_wvlen_ratio = 0.2;   % vary 20% from the basic L0;
[ts_pair, basic_parms] = create_simple_signals_with_rednoise(ncount, x, delta_phase, delta_wvlen_ratio);
save(['red_noise_signal_n' num2str(ncount) '_v2.mat'],'ts_pair','basic_parms');


load('create_demo/red_noise_signal_n100_v2.mat');
ATOMIC_platform = 'sine wave w/ red noise';
xres = 0.2;
hw = waitbar(0, 'Running wavelet coherence...');
for i = 10:ncount
    
    % perform wavelet coherence computation:
    testobj = ATOMIC_Wavelet(ts_pair(i).ts1, ts_pair(i).ts2, x, time, ATOMIC_platform, xres);
    
    [hfig, wtc_stat(i), outCOI_stat(i), inCOI_stat(i)]=testobj.wavelet_coherence_toolbox;

    if mod(i, floor(ncount/10))<0.001
        waitbar(i/ncount, hw);
        
        xc_savefig(hfig, figsvdir,['sine wave with random red noise instance' num2str(i,'%3.3i') '.jpg'],[0 0 8 6]);
        
        labelstr.ts1 = {'fake SST', 'w/ red noise'};
        labelstr.ts2 = {'sine wave', 'wind speed'};
        testobj.plot_data(labelstr);
        xc_savefig(gcf, figsvdir,['sample_signal_instance' num2str(i,'%3.3i') '.jpg'],[0 0 8 6]);
        
    end

    
end

save(['./create_demo/wavelet_coherence_for_100_tspair_v2.mat'],'wtc_stat','outCOI_stat','inCOI_stat');


%% now plot the results. (try establishing another object so that I don't need to repeat this code too many times.
demo_obj = ATOMIC_Wavelet(ts_pair(i).ts1, ts_pair(i).ts2, x, time, ATOMIC_platform, xres);
ds_out.wvlet_coh = wtc_stat;
ds_out.outCOI_metrics = outCOI_stat;
ds_out.inCOI_metrics = inCOI_stat;
demo_obj.wtcstat = ds_out;

inCOI_flag = false;
hfig0= demo_obj.plot_scale_and_phase(99,inCOI_flag);
hfig0.Position = hfig.Position;
rlim([0 45])
xc_savefig(gcf, figsvdir, ['phase_wavelength_scatter_plots_100_rednoise_sinewaves_pair_v2_GoodOnly.jpg'], [0 0 10,8]);


% now, I can use kcluster to find out the groups I visually identified, and
% get the averaged wavelenth from the group, as well as the standard
% deviation as a guess of the true value;
Xtest = [[outCOI_stat.ave_phase]', [outCOI_stat.ave_wvlen]'];
[kidx, C]= kmeans(Xtest,3);
figure;
for i= 1:3
polarscatter(Xtest(kidx==i,1), Xtest(kidx==i,2),'.');
%'.b');
hold on
end
th = [0:0.05*pi:2*pi];
polarplot(th,P90s.QC*ones(size(th)),'--k','linewidth',1.2);
%polarscatter(Xtest(kidx==2,1), Xtest(kidx==2,2),'.r');
title('Cluster Groups by k-means clustering');
set(gca,'fontsize',14);

% use visual cue and wave length for separation:
true_guess_crit = [outCOI_stat.ave_wvlen]>20;
true_metric_g = Xtest(true_guess_crit,:);
true_metric_gm = mean(true_metric_g,1);
true_metric_gstv =std(true_metric_g,0,1); % w default is 0, averaged by (N-1) observations.

hold on;
polarplot(true_metric_gm(1), true_metric_gm(2), '*m','markersize',12,'linewidth',2);
ph_range_g = true_metric_gm(1)+[-true_metric_gstv(1):0.01*pi:true_metric_gstv(1)];
polarplot(ph_range_g, true_metric_gm(2).*ones(size(ph_range_g)),':m','linewidth',2);
% plot oscillation range in wavelength:
polarplot([true_metric_gm(1), true_metric_gm(1)], ...
    true_metric_gm(2)+[-true_metric_gstv(2),true_metric_gstv(2)],':m','linewidth',2);
xc_savefig(gcf, figsvdir, ['phase_wavelength_scatter_plots_100_rednoise_sinewaves_pair_v2_GoodOnly_guessedAnswer.jpg'], [0 0 10,8]);


inCOI_flag = true;
hfig= demo_obj.plot_scale_and_phase(99,inCOI_flag);
rlim([0 45])
xc_savefig(hfig, figsvdir, ['phase_wavelength_scatter_plots_100_rednoise_sinewaves_pair_v2_withEdgeEffects.jpg'], [0 0 10,8]);
% mark the true value on the chart;
hold on;
% negative because ts2 is leading ts1.
polarplot(-basic_parms.phase, basic_parms.wvlen, 'pm','markersize',15,'linewidth',2);
% plot oscillation range in phase:
ph_range = -[basic_parms.phase-delta_phase:0.01*pi:basic_parms.phase+delta_phase];
polarplot(ph_range, basic_parms.wvlen.*ones(size(ph_range)),'-m','linewidth',2);
% plot oscillation range in wavelength:
polarplot(-[basic_parms.phase, basic_parms.phase], ...
    basic_parms.wvlen.*[1-delta_wvlen_ratio, 1+delta_wvlen_ratio],'-m','linewidth',2);

xc_savefig(gcf, figsvdir, ['phase_wavelength_scatter_plots_100_rednoise_sinewaves_pair_v2_withEdgeEffects_addedTrueAnswer.jpg'], [0 0 10,8]);

% add the guessed value on the same plot:
hold on
polarplot(true_metric_gm(1), true_metric_gm(2), '*k','markersize',12,'linewidth',2);
ph_range_g = true_metric_gm(1)+[-true_metric_gstv(1):0.01*pi:true_metric_gstv(1)];
polarplot(ph_range_g, true_metric_gm(2).*ones(size(ph_range_g)),':k','linewidth',2);
% plot oscillation range in wavelength:
polarplot([true_metric_gm(1), true_metric_gm(1)], ...
    true_metric_gm(2)+[-true_metric_gstv(2),true_metric_gstv(2)],':k','linewidth',2);
xc_savefig(gcf, figsvdir, ['phase_wavelength_scatter_plots_100_rednoise_sinewaves_pair_v2_all.jpg'], [0 0 10,8]);


inCOI_flag = true;
[hfig, P90s]= demo_obj.plot_scale_and_phase(99,inCOI_flag);
