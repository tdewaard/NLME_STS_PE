function [beta,PSI,stats,b,t_elapsed] = nlme_param_est(model,time,data,groups,graphics)
%NLME_PARAM_EST Summary of this function goes here
%   Detailed explanation goes here
if graphics
    outfcn = @nlmefitoutputfcn;
else
    outfcn = [];
end
disp(strcat("Beginning NLME parameter estimation with model ", model.name,"..."))
statsetopts=statset('Display', 'iter', 'FunValCheck', 'on','OutputFcn', outfcn);
tic
try
    [beta,PSI,stats,b] = nlmefit(time,data,groups, [], model.nlmehandle, log(model.p0),...
                              'Options', statsetopts,...
                              'RefineBeta0','on',...            
                              'ErrorModel', 'constant',...                          %error model y = f + a*error
                              'ParamTransform', zeros(1,length(model.p0)),...       %param transform, ones=log conversion, zeros=no conversion NOT NEEDED
                              'CovPattern', ones(length(model.p0)));
catch
    switch model.name
        case 'RIBO'
            model.nlmehandle=@(PHI,t)nlme_model_ribo(PHI,t,true,true);
        case 'ENZDEGRIBO'
            model.nlmehandle=@(PHI,t)nlme_model_enzdegribo(PHI,t,true,true);
   end
   [beta,PSI,stats,b] = nlmefit(time,data,groups, [], model.nlmehandle, log(model.p0),...
                              'Options', statsetopts,...
                              'RefineBeta0','on',...            
                              'ErrorModel', 'constant',...                          %error model y = f + a*error
                              'ParamTransform', zeros(1,length(model.p0)),...       %param transform, ones=log conversion, zeros=no conversion NOT NEEDED
                              'CovPattern', ones(length(model.p0)));
end
t_elapsed=toc;
disp("NLME Parameter estimation completed.")
end

