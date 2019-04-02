function [pest,resnorm,resid,exitflag,output,lambda,Jacobian,covarpar,t_elapsed] = param_est(model,p0,texp,fitdata,genplots,totaldata)
%PARAM_EST Summary of this function goes here
%   Detailed explanation goes here
lb = zeros(size(model.p0));
ub = [];
options = optimoptions('lsqnonlin',...      %default options
                       'Display','final',...
                       'FunctionTolerance', 1e-6,...
                       'MaxFunctionEvaluations', 800,...
                       'StepTolerance', 1e-6,...
                       'OptimalityTolerance',1e-6);

switch model.name
    case 'NORM'
        errorfun=@transfection_norm_error;
    case 'RIBO'
        errorfun=@transfection_ribo_error;
end

disp(strcat("Beginning parameter estimation with ",model.name," model..."))
tic
[pest,resnorm,resid,exitflag,output,lambda,Jacobian] = lsqnonlin(errorfun, p0, lb, ub, options, texp, fitdata);
t_elapsed=toc;

disp(strcat("Parameter estimation completed after ",num2str(t_elapsed)," seconds."))

Jacobian = full(Jacobian);
covarpar = resnorm*inv(Jacobian'*Jacobian)/length(fitdata);

%plotting
if genplots
    ci = nlparci(pest,resid,'jacobian',Jacobian);

    %Beta coefficient plot with confidence intervals
    %Dot-and-Whisker Plots of Regression Results
    figure('Name','parameter estimates for population average') 

    errorbar(log10(pest),[1:length(pest)],log10(pest')-log10(ci(:,1)), 'horizontal','ob', 'LineWidth',2); hold on

    ylim([0 length(model.pnames)+1])
    yticks(1:length(model.pnames))
    yticklabels(model.pnames)
    ax = gca; ax.YDir = 'reverse'; ax.FontSize = 16;
    xlabel('Log10')
    title(strcat('Parameter estimates of model ',model.name,' with 95% confidence interval'))

    %plot calibrated model and data
    [h1,h2]=make_data_plot(totaldata,texp,true);
    hold on
    [t,x,y]=simulate_model(model,texp,[],pest,false,2);
    h3=plot(t,y,'b', 'LineWidth',3);
    hold off
    legend([h1, h2(1), h3],'average data', '+/- standard deviation of the data','fitted model')
end
end

