function [ts_pair, basic_parm] = create_simple_signals_with_rednoise(ncount,x,delta_phase, delta_wvlen_ratio)

% ncount: number of sample sets
% x: input distance or time axis

% blindly generate the base phase and wavelength for the time series;
% ph0 = rand(1);  % require it to be within 0 and 2*pi;
% L0 = rand(1);   %require it to be within 10 to 30km

[ph0, L0] = generate_basic_parameter;
LX = max(x) - min(x);
delta_wvlen = delta_wvlen_ratio*L0;

hw = waitbar(0,'running..');

for i = 1:ncount
    % generate noise (oscilation in phase and wvlen within the given
    % range);
    if mod(i,2)==0
        sign_here = 1;
    else
        sign_here = -1;
    end
    theta = ph0 + sign_here*rand(1)* delta_phase;
    wvlen = L0 + sign_here*rand(1)* delta_wvlen;
    
    x0_crit = false;
    while ~x0_crit
        rtmp = rand(1);
        x0_crit = (rtmp>0.2 & rtmp<0.8);
    end
    x0 = rtmp*LX;
    
    k1 = 2*pi/wvlen;
    
    % generate red noise for the two series:
    rrn = generate_red_noise;
    
    ts1_tmp = sin(k1.* (x - x0));
    ts1_tmp(x<(x0+0.25*wvlen-0.75*wvlen))=0;
    ts1_tmp(x>(x0+0.25*wvlen+0.75*wvlen))=0;
    
    ts_pair(i).ts1 = ts1_tmp +  rrn(:,1)';
    ts_pair(i).ts2 = sin((k1.*(x - x0) + theta)) + rrn(:,2)';
    ts_pair(i).ts3 = sin((k1.*(x - x0) + theta)) + sin((0.5*k1.*(x - x0) + theta+pi/2)) + rrn(:,2)';
    
    
    if mod(i, floor(ncount/10))<0.01
        waitbar(i/ncount, hw);
        
        figure(10);
        plot(x, ts_pair(i).ts1,'-b');
        hold on
        plot(x, ts_pair(i).ts2,'-r');
        plot(x, ts_pair(i).ts3,'-m');
        hold off
        pause
    end
end

basic_parm.phase = ph0;
basic_parm.wvlen = L0;



%%%%%%%%%%%%%%%%%%% Nested Functions %%%%%%%%%%%%%%%%%%%
    function [ph0, L0] = generate_basic_parameter()

        % initialize:
        ph_crit = false;
        L0_crit = false;
       
        while ~ph_crit
            ph0_tmp = rand(1)*2*pi;
            ph_crit = ph0_tmp>=0 & ph0_tmp<=2*pi;
        end
        
        while ~L0_crit
            L0_tmp = rand(1)*100;
            L0_crit= L0_tmp >10 & L0_tmp <35;
        end
        
        ph0 = ph0_tmp;
        L0 = L0_tmp;
        
        disp('finished generating basic parameters.');
        
    end

    % try with the matlab object dsp.ColoredNoise
    function redn = generate_red_noise()
       cn =  dsp.ColoredNoise('brown','SamplesPerFrame',length(x), 'NumChannels',2, ...
            'BoundedOutput',false);
       redn = cn();
      
       redn(:,1) = redn(:,1)./max(abs(redn(:,1)));
       redn(:,2) = redn(:,2)./max(abs(redn(:,2)))*2;
       
       figure(11);
       subplot(2,1,1)
       plot(redn(:,1)); title('Channel 1'); axis tight;
       subplot(2,1,2)
       plot(redn(:,2)); title('Channel 2'); axis tight;
    end





end