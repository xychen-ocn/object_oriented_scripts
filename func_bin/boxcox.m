% boxcox transformation implementation:
function [transdat, lambda_optimal, maxllf] = boxcox(data)

% maximizie the log-likelihood ==> minimize the negative log likelihood;
%data = x; 
if ~iscolumn(data)
    data = data';
end


lambda_guess = linspace(-5,100);
llf=[];
for ii = 1:length(lambda_guess)
    lambda = lambda_guess(ii);    
    y = boxcox_transfunc(lambda, data);
    llf(ii) = boxcox_llf(lambda, data,y);
end
figure(11)
plot(lambda_guess, llf);
set(gca,'yscale','log');
pause(0.5)

[maxllf, mid] = max(llf);
lambda_optimal = lambda_guess(mid);
transdat = boxcox_transfunc(lambda_optimal, data);


% options = optimset('PlotFcns',@optimplotfval);
% 
% lambda = fminsearch(@myfun, 0); %, options);


% function nLL = myfun(lambda)
% disp(num2str(lambda));
% if lambda ==0
%     transdat = log(data);
% else
%     transdat = (data.^lambda -1)./lambda;
% end
% 
% % find the negative log likelyhood
% % if mean(transdat)<0
% %     transdat = transdat+mean(transdat);
% % end
% pd = fitdist(transdat,'normal');
% nLL = (negloglik(pd));
% 
% end


    function llf = boxcox_llf(lambda, x, y)
        % x: input, y: transformed x;
        N = length(x);
        llf = (lambda - 1).* sum(log(x)) - ...
        N/2 * log(sum((y - mean(y)).^2)./N);
    end

    function y = boxcox_transfunc(lambda, data)
        if lambda ==0
            y = log(data);
        else
            y = (data.^lambda -1)./lambda;
        end
    end

end


