%%
%Script for processing results / making figures
%
%
%

clear; close all

addpath('functions',genpath('Results'))

MODELS = init_models();                         %MODELS(1) = NORM, MODELS(2) = RIBO

dates = ["20160427","20160513", "20160613"];
[texp,fluor] = extract_data(dates(2),'stable',true);

results.STS.NORM = dir(fullfile('Results\LSQNONLIN\NORM','*.mat'));
results.STS.RIBO = dir(fullfile('Results\LSQNONLIN\RIBO','*.mat'));
results.NLME.NORM.Real = dir(fullfile('Results\NLMEFIT\NORM\RealData','*.mat'));
results.NLME.NORM.Synth = dir(fullfile('Results\NLMEFIT\NORM\SynthData','*.mat'));
results.NLME.RIBO = dir(fullfile('Results\NLMEFIT\RIBO','*.mat'));

%% Raw Data Plot

make_data_plot(fluor,texp,true);

%% All Models LSQNONLIN NORM

load(results.STS.NORM(1).name)
figure; hold on
title("\bf{All Models LSQNONLIN NORM}")
xlim([0,texp(end)])
for c=1:length(Out)
    p=Out{c}.Theta;
    simulate_model(MODELS(1),texp,[],p,true,[0 0 0 0.1],1);
end
hold off

%% All Models LSQNONLIN RIBO

load(results.STS.RIBO(1).name)
figure; hold on
title("\bf{All Models LSQNONLIN RIBO}")
xlim([0,texp(end)])
for c=1:length(Out)
    p=Out{c}.Theta;
    simulate_model(MODELS(2),texp,[],p,true,[0 0 0 0.1],1);
end
hold off

%% 10 Random Model/Data LSQNONLIN NORM

load(results.STS.NORM(1).name)
r = randi([1,length(fluor)],10,1)';
figure; hold on
title("\bf{10 Random Model/Data LSQNONLIN NORM}")
xlim([0,texp(end)])
for k = r
    p=Out{k}.Theta;
    simulate_model(MODELS(1),texp,[],p,true,'k',2);
    plot(texp,fluor(:,k),'--r','Linewidth',2);
end

%% 10 Random Model/Data LSQNONLIN RIBO

load(results.STS.RIBO(1).name)
r = randi([1,length(fluor)],10,1)';
figure; hold on
title("\bf{10 Random Model/Data LSQNONLIN RIBO}")
xlim([0,texp(end)])
for k = r
    p=Out{k}.Theta;
    simulate_model(MODELS(2),texp,[],p,true,'k',2);
    plot(texp,fluor(:,k),'--r','Linewidth',2);
end

%% 10 Random Model/Data NLMEFIT NORM RealData

n_groups=length(fluor);

for i=1:1:length(results.NLME.NORM.Real)
    n=results.NLME.NORM.Real(i).name;
    if contains(n,strcat(num2str(n_groups),'G'))
        Out = load(n);
        disp("NLME dataset loaded!")
    end
end

PHI_norm = repmat(Out.Beta,1,n_groups)+Out.b;

figure; hold on
xlim([0,texp(end)])
title("10 Random Model/Data NLMEFIT NORM")

r = randi([1,n_groups],10,1)';
for i=r
    y = MODELS(1).nlmehandle(PHI_norm(:,i)',texp);
    mod=plot(texp,y,'k','Linewidth',2);
    dat=plot(texp,fluor(:,i),'r--','Linewidth',2);
end
legend([mod,dat],{"NORM Model Trace","Experimental Data Trace"}, 'location','southeast')
ylabel("Fluorescence (a.u.)")
xlabel("Time (h)")
set(findall(gcf,'-property','FontSize'),'FontSize',13)

%% LSQNONLIN Elapsed time distribution 

n = load(results.STS.NORM(1).name);
r = load(results.STS.RIBO(1).name);
t_NORM=zeros(length(fluor),1);
t_RIBO=zeros(length(fluor),1);
for i=1:length(fluor)
    t_NORM(i) = n.Out{i}.ElapsedTime;
    t_RIBO(i) = r.Out{i}.ElapsedTime;
end
figure;
subplot(2,1,1)
histogram(t_NORM,"BinWidth",0.01,"EdgeAlpha",0.2,'Facecolor',[1 0 0.1])
title(strcat("\bf{Computation time NORM model, \mu = ",num2str(round(mean(t_NORM),3)),", ",num2str(round(sum(t_NORM),3))," seconds total}"))
xlabel("Time (s)")

subplot(2,1,2)
histogram(t_RIBO,"BinWidth",0.07,"EdgeAlpha",0.2,'Facecolor',[0.1 1 0.1])
title(strcat("\bf{Computation time RIBO model, \mu = ",num2str(round(mean(t_RIBO),3)),", ",num2str(round(sum(t_RIBO),3))," seconds total}"))
xlabel("Time (s)")
set(findall(gcf,'-property','FontSize'),'FontSize',13)

%t-test, to verify RIBO has significantly higher computation times.
[h,p_value,ci,stats]=ttest(t_RIBO,t_NORM,'Tail','right')

%% LSQNONLIN Parameter distribution 
r=load(results.STS.RIBO(1).name);
n=load(results.STS.NORM(1).name);

Theta_norm=zeros(length(fluor),length(MODELS(1).p0));
Theta_ribo=zeros(length(fluor),length(MODELS(2).p0));
for k=1:length(n.Out)
    Theta_norm(k,:)=n.Out{k}.Theta;
    Theta_ribo(k,:)=r.Out{k}.Theta;
end

%norm distribution
spn=321;
figure;
sgtitle('\bf{Parameter distributions NORM model}')
for i=1:length(MODELS(1).p0)
    subplot(spn)
    m=round(mean(Theta_norm(:,i)),3);
    histogram(Theta_norm(:,i),"BinWidth",0.05*m,"EdgeAlpha",0.2,'Facecolor',[1 0 0.1]);hold on
    Mean_line=line([m, m], get(gca, 'ylim'),'LineStyle', '-.', 'LineWidth', 2, 'Color', 'b');
    xlabel(MODELS(1).pnames{i})
    title(strcat("\bf{\mu = ",num2str(m),"}"))
    spn=spn+1;
end
hold off
legend([Mean_line],{"Parameter mean"})
set(findall(gcf,'-property','FontSize'),'FontSize',13)

%ribo distribution
spn=421;
figure;
sgtitle('\bf{Parameter distributions RIBO model}')
for i=1:length(MODELS(2).p0)
    subplot(spn)
    m=round(mean(Theta_ribo(:,i)),3);
    histogram(Theta_ribo(:,i),"BinWidth",0.05*m,"EdgeAlpha",0.2,'Facecolor',[0.1 1 0.1]);hold on
    Mean_line=line([m, m], get(gca, 'ylim'),'LineStyle', '-.', 'LineWidth', 2, 'Color', 'b');
    xlabel(MODELS(2).pnames{i})
    title(strcat("\bf{\mu = ",num2str(m),"}"))
    spn=spn+1;
end
hold off
legend([Mean_line],{"Parameter mean"})
set(findall(gcf,'-property','FontSize'),'FontSize',13)

%% NLMEFIT Parameter distribution

%extract NLME data
n_groups=length(fluor);
for i=1:1:length(results.NLME.NORM.Real)
    n=results.NLME.NORM.Real(i).name;
    if contains(n,strcat(num2str(n_groups),'G'))
        Out = load(n);
        disp("NLME dataset loaded!")
    end
end

PHI_norm = exp(repmat(Out.Beta,1,n_groups)+Out.b);   %mixed effects
Beta_norm = exp(Out.Beta);                           %fixed effects

figure;
sgtitle('\bf{Mixed Effect Parameter distributions NORM model}')
hold on
for i=1:length(MODELS(1).p0)
    subplot(3,2,i)
    m=round(mean(PHI_norm(i,:)),3);
    histogram(PHI_norm(i,:),"BinWidth",0.05*m,"EdgeAlpha",0.2,'Facecolor',[1 0 0.1]);
    B_line=line([Beta_norm(i), Beta_norm(i)], get(gca, 'ylim'),'LineStyle', '--', 'LineWidth', 2, 'Color', 'k');
    Mean_line=line([m, m], get(gca, 'ylim'),'LineStyle', '-.', 'LineWidth', 2, 'Color', 'b');
    xlabel(MODELS(1).pnames{i})
    title(strcat("\bf{\mu = ",num2str(m),", \beta = ",num2str(Beta_norm(i)),"}"))
end
legend([B_line,Mean_line],{"Fixed effects \beta", "Parameter mean"})
hold off
set(findall(gcf,'-property','FontSize'),'FontSize',13)

%% NLMEFIT ElapsedTime + RMSE Distribution

Out=cell(length(results.NLME.NORM.Real),1);
t=zeros(length(results.NLME.NORM.Real),1);
g=zeros(length(results.NLME.NORM.Real),1);

n_groups=length(fluor);

for i=1:1:length(results.NLME.NORM.Real)
    n=results.NLME.NORM.Real(i).name;
    if contains(n,strcat(num2str(n_groups),'G'))
        O = load(n);
        disp("NLME dataset loaded!")
    end
end

PHI_norm = repmat(O.Beta,1,n_groups)+O.b;

%calculate all RMSE's
rmse=zeros(1,length(fluor));
for i=1:n_groups
    y = MODELS(1).nlmehandle(PHI_norm(:,i)',texp);
    rmse(i)=calc_rmse(fluor(:,i),texp,y);
end

for n=1:n_groups
    for i=1:1:length(results.NLME.NORM.Real)
        fname=results.NLME.NORM.Real(i).name;
        if contains(fname,strcat(num2str(n),'G'))
            Out{i} = load(fname);
            t(i)=Out{i}.t_elapsed;
            g(i)=n;
            disp(strcat("NLME dataset with ",num2str(n)," groups loaded!"))
        end
    end
end

%fitfun=@(a,x)a(1)*log(x+a(2))+a(3);
%[a,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(fitfun,[1 1 1],g,t);
%RMSE=sqrt(resnorm)/length(g);
[h,p_value,ci,stats]=ttest(rmse,RMSE_norm','Tail','left')


subplot(2,1,1);
    loglog(g,t,'o','MarkerSize',6,'Linewidth',3); hold on;
    %loglog(linspace(2,1000),fitfun(a,linspace(2,1000)),'k-',"LineWidth",2)
    title("\bf{Computation times NLME model NORM}")
    xlabel("Number of individuals")
    ylabel("Elapsed time during NLME (s)")
    grid on
subplot(2,1,2);
    m=mean(rmse);
    histogram(rmse,"BinWidth",0.005,"EdgeAlpha",0.2,'Facecolor',[1 0 0.1])
    title(strcat("\bf{NLME Result NORM RMSE Distribution, \mu = ",num2str(m),"}"))
    xlabel("RMSE (a.u.)")
set(findall(gcf,'-property','FontSize'),'FontSize',13)

%% STS MeanData NORM + RIBO

n=load(results.STS.NORM(2).name);
r=load(results.STS.RIBO(2).name);

%params
theta_norm=n.s.Theta;
theta_ribo=r.s.Theta;

%param covariances
resid_norm=n.s.Residuals;
resid_ribo=r.s.Residuals;

jac_norm = n.s.Jacobian;
jac_ribo = r.s.Jacobian;

ci_norm = nlparci(theta_norm,resid_norm,'jacobian',jac_norm);
ci_ribo = nlparci(theta_ribo,resid_ribo,'jacobian',jac_ribo);

%performance measure elapsed time and RMSE
t_elapsed=[n.s.ElapsedTime r.s.ElapsedTime];
rss=[n.s.Resnorm r.s.Resnorm];
RMSE=sqrt(rss/length(texp));

%plotting
subplot(2,2,1)
    xlim([0,texp(end)])
    title("\bf{Data With Mean and NORM Model Fit}")
    hold on
        plot(texp,fluor,'Color',[0 0 0 0.1])
        simulate_model(MODELS(1),texp,[],theta_norm,true,'r',2);
        plot(texp,mean(fluor')','c--',"Linewidth",2)
    hold off
    text(10,4,strcat("\it{RMSE = ",num2str(round(RMSE(1),3)),"}"))
    text(10,3,strcat("\it{t_{elapsed} = ",num2str(round(t_elapsed(1),3))," s}"))
subplot(2,2,2)
    errorbar(theta_norm,[1:length(theta_norm)],theta_norm'-ci_norm(:,1), 'horizontal','ob', 'LineWidth',2);
    ylim([0 length(MODELS(1).pnames)+1])
    yticks(1:length(MODELS(1).pnames))
    yticklabels(MODELS(1).pnames)
    ax = gca; ax.YDir = 'reverse';
    xlabel('ln(\theta_{mean})')
    title("\bf{Parameter Estimates NORM model 95% confidence interval}")
subplot(2,2,3)
    xlim([0,texp(end)])
    title("\bf{Data With Mean and RIBO Model Fit}")
    hold on
        plot(texp,fluor,'Color',[0 0 0 0.1])
        simulate_model(MODELS(2),texp,[],theta_ribo,true,'r',2);
        plot(texp,mean(fluor')','c--',"Linewidth",2)
    hold off
    text(10,4,strcat("\it{RMSE = ",num2str(round(RMSE(2),3)),"}"))
    text(10,3,strcat("\it{t_{elapsed} = ",num2str(round(t_elapsed(2),3))," s}"))
subplot(2,2,4)
    errorbar(theta_ribo,[1:length(theta_ribo)],theta_ribo'-ci_ribo(:,1), 'horizontal','ob', 'LineWidth',2);
    ylim([0 length(MODELS(2).pnames)+1])
    yticks(1:length(MODELS(2).pnames))
    yticklabels(MODELS(2).pnames)
    ax = gca; ax.YDir = 'reverse';
    xlabel('ln(\theta_{mean})')
    title("\bf{Parameter Estimates RIBO model 95% confidence interval}")
set(findall(gcf,'-property','FontSize'),'FontSize',13)

%% STS NORM + RIBO
n=load(results.STS.NORM(1).name);
r=load(results.STS.RIBO(1).name);

RMSE_norm=zeros(length(fluor),1);
RMSE_ribo=zeros(length(fluor),1);
for i=1:length(fluor)
    rss_norm=n.Out{i}.Resnorm;
    rss_ribo=r.Out{i}.Resnorm;
    RMSE_norm(i)=sqrt(rss_norm./length(texp));
    RMSE_ribo(i)=sqrt(rss_ribo./length(texp));
end
[h,p_value,ci,stats]=ttest(RMSE_ribo,RMSE_norm,'Tail','left')

subplot(2,3,1)
    hold on
        title("\bf{All Models LSQNONLIN NORM}")
        xlim([0,texp(end)])
        for c=1:length(n.Out)
            p=n.Out{c}.Theta;
            simulate_model(MODELS(1),texp,[],p,true,[0 0 0 0.1],1);
        end
    hold off
    
subplot(2,3,2)
    m=round(mean(RMSE_norm),3);
    histogram(RMSE_norm,"BinWidth",0.005,"EdgeAlpha",0.2,'Facecolor',[1 0 0.1])
    title(strcat("\bf{RMSE Distribution NORM, \mu = ",num2str(m),"}"))
    xlabel("RMSE (a.u.)")
    
subplot(2,3,3)
    hold on
    title("\bf{10 Random Model/Data LSQNONLIN NORM}")
    xlim([0,texp(end)])
    for k = randi([1,length(fluor)],10,1)'
        p=n.Out{k}.Theta;
        simulate_model(MODELS(1),texp,[],p,true,'k',2);
        plot(texp,fluor(:,k),'--r','Linewidth',2);
    end
    hold off
    
subplot(2,3,4)
    hold on
        title("\bf{All Models LSQNONLIN RIBO}")
        xlim([0,texp(end)])
        for c=1:length(r.Out)
            p=r.Out{c}.Theta;
            simulate_model(MODELS(2),texp,[],p,true,[0 0 0 0.1],1);
        end
    hold off
    
subplot(2,3,5)
    m=round(mean(RMSE_ribo),3);
    histogram(RMSE_ribo,"BinWidth",0.005,"EdgeAlpha",0.2,'Facecolor',[0.1 1 0.1])
    title(strcat("\bf{RMSE Distribution RIBO, \mu = ",num2str(m),"}"))
    xlabel("RMSE (a.u.)")
    
subplot(2,3,6)
    hold on
    title("\bf{10 Random Model/Data LSQNONLIN RIBO}")
    xlim([0,texp(end)])
    for k = randi([1,length(fluor)],10,1)'
        p=r.Out{k}.Theta;
        simulate_model(MODELS(2),texp,[],p,true,'k',2);
        plot(texp,fluor(:,k),'--r','Linewidth',2);
    end
    hold off
set(findall(gcf,'-property','FontSize'),'FontSize',13)

%% NLMEFIT NORM Mean Characteristics Comparisons

%loading in necessary data for comparison
norm_all=load(results.STS.NORM(1).name);
norm_mean=load(results.STS.NORM(2).name);
theta_norm_mean=norm_mean.s.Theta;

%extract all the param estimates from STS
p=zeros(length(norm_all.Out),length(MODELS(1).p0));
for c=1:length(norm_all.Out)
    p(c,:)=norm_all.Out{c}.Theta;
end
mean_params_STS=mean(p,1);   %mean parameters

%extract NLME data
n_groups=length(fluor);
for i=1:1:length(results.NLME.NORM.Real)
    n=results.NLME.NORM.Real(i).name;
    if contains(n,strcat(num2str(n_groups),'G'))
        Out = load(n);
        disp("NLME dataset loaded!")
    end
end

PHI_norm = repmat(Out.Beta,1,n_groups)+Out.b;   %mixed effects
Beta_norm = repmat(Out.Beta,1,n_groups);        %fixed effects

%plotting
figure; hold on
xlim([0,texp(end)])
title(strcat("NORM model simulation NLME with ",num2str(n_groups)," individuals"))
for i=1:n_groups
    y = MODELS(1).nlmehandle(PHI_norm(:,i)',texp);
    y_fixed = MODELS(1).nlmehandle(Beta_norm(:,i)',texp);
    h1=plot(texp,y,'Color',[0 0 0 0.1],'Linewidth',1);
    h2=plot(texp,y_fixed,'r','linewidth',3);
end
[~,~,y_mean,h4]=simulate_model(MODELS(1),texp,[],theta_norm_mean,true,[0.5 1 0.5],3);
[~,~,y_mean_param_STS,h5]=simulate_model(MODELS(1),texp,[],mean_params_STS,true,[1 0.5 1],3);
h3=plot(texp, mean(fluor')','c--','linewidth',3);

%calc necessary rmse's
rmse=round([calc_rmse(mean(fluor')',texp,y_fixed), calc_rmse(mean(fluor')',texp,y_mean),calc_rmse(mean(fluor')',texp,y_mean_param_STS)],3);

legend([h3,h1,h2,h4,h5],{"Mean data", ...
                        "NORM model mixed effects simulations",...
                        strcat("NORM model fixed effects only, \it{RMSE} = ",num2str(rmse(1))),...
                        strcat("NORM model calibrated on Mean, \it{RMSE} = ",num2str(rmse(2))),...
                        strcat("NORM model w/ mean STS params, \it{RMSE} = ",num2str(rmse(3)))},...
                        'Location','southeast')
hold off
set(findall(gcf,'-property','FontSize'),'FontSize',13)




