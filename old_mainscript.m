


%% Loading in data
clear; close all
addpath('functions')

dates = ["20160427","20160513", "20160613"];

[texp,fluor] = extract_data(dates(2),'stable',true);

%% Plotting data

make_data_plot(fluor,texp,true);

%% initialise models
MODELS = init_models();

%% plotting of uncalibrated models
figure; hold on;
simulate_model(NORM,texp,[],NORM.p0,true,'r',3);
simulate_model(RIBO,texp,[],RIBO.p0,true,'b',3);
legend('NORM model','RIBO model')
title("\bf{Uncalibrated model simulation}")
xlim([0,texp(end)])
hold off

%% LSQNONLIN Parameter estimation with mean data.

%output automatically saved to file
for m=1:1:length(MODELS)
    [pest,resnorm,resid,exitflag,output,lambda,Jacobian,covarpar,t_elapsed] = param_est(MODELS(m),MODELS(m).p0,texp,mean(fluor')',false, fluor);
    s = struct('Theta',pest,'Resnorm',resnorm,'Residuals',resid,'Exitflag',exitflag,'Output',output,'Lambda',lambda,'Jacobian',Jacobian,'CovMatrix',covarpar,'ElapsedTime',t_elapsed);
    filename=strcat("MeanDataResult_LSQNONLIN_",MODELS(m).name,".mat");
    save(filename,'s')
end

%% LSQNONLIN Parameter estimation with all data, all models

for m=1:1:length(MODELS)
    Out=cell(length(fluor),1);
    h=waitbar(0,strcat("Please wait, calibrating ",MODELS(m).name," models..."));
    for k = 1:length(fluor)
        disp(strcat("Calibrating model ",MODELS(m).name,", individual #",num2str(k)," out of ",num2str(length(fluor))))
        [pest,resnorm,resid,exitflag,output,lambda,Jacobian,covarpar,t_elapsed] = param_est(MODELS(m),MODELS(m).p0,texp,fluor(:,k),false,fluor);
        Out{k,:} = struct('Theta',pest,'Resnorm',resnorm,'Residuals',resid,'Exitflag',exitflag,'Output',output,'Lambda',lambda,'Jacobian',Jacobian,'CovMatrix',covarpar,'ElapsedTime',t_elapsed);
        disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ITERATION DONE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        h=waitbar(k / length(fluor)); 
    end
    close(h)
    disp(strcat("Calibrating with model ",MODELS(m).name," done, saving output..."))
    save(strcat('AllDataResult_LSQNONLIN_',MODELS(m).name,'.mat'),'Out')
    disp("==========================================================================================")
end

%% Some plotting

load('Pest') 

spn=421;
for i=1:size(Pest,2)
    subplot(spn)
    histfit( log10(Pest(:,i)),10,'kernel'); hold on
    xlabel(pnames{i})
    spn=spn+1;
end

%plot all models
figure; hold on;
for k = 1:length(Pest)
    disp(strcat("Simulating model of individual #",num2str(k)))
    [t,x,y] = simulate_model(RIBO,texp,[],Pest(k,:),true,[0 0 0 0.1],1); 
end
xlim([0,texp(end)])
ylim([0,max(fluor(:))])
title('all models')
hold off

%plot some models and corresponding data
figure
r = randi([1,length(fluor)],10,1)';
for k = r
    [t,x,y] = simulate_model(RIBO,texp,[],Pest(k,:),false,'k',1);
    plot(texp,fluor(:,k),'--r'); hold on
    plot(t,y,'b');
end
xlabel('Time (h)')
ylabel('Fluorescence Intensity (a.u.)')
legend('Data','Model')
title('\bf{10 Random Exxamples}')


%% ========================== NLME ==========================

save("Results\test.mat",'b')
%% Generate random synthetic data

%mixed-effect on synth data

for G=6:10
    m=RIBO;
    n_groups=G;
    sigma_xi=0.05;
    
    [sTIME, sGROUPS, sDATA]=generate_synth_data(m,sigma_xi,1,false,n_groups,false);
    %NLMEfit on synthetic data
    [Beta,Psi,stats,b,t_elapsed]=nlme_param_est(m,sTIME,sDATA,sGROUPS,false)
    filename=strcat("Results\NLMEFIT\",m.name,"\SynthData",num2str(n_groups),"G_xi=",num2str(sigma_xi),"Result_NLMEFIT_",m.name,".mat");
    if ~isfile(filename)
        save(filename,'Beta','Psi','stats','b','t_elapsed')
    else
        save(strcat("NLMEtestXI=",num2str(sigma_xi),"G=",num2str(n_groups),".mat"),'Beta','Psi','stats','b','t_elapsed')
    end
end
%%
%convert real data to 'grouped data' format for NLMEfit
TIME=repmat(texp,length(fluor),1)';
TIME=TIME(:);
GROUPS=repmat(1:length(fluor),length(texp),1);
GROUPS=GROUPS(:);
DATA=fluor(:);

%% Post processing NLME results
norm_results = load("NLMEresults_norm.mat");
ribo_results = load("NLMEresults_ribo.mat");

%combine fixed effect with random effect in PHI
PHI_norm = repmat(norm_results.phi,1,6)+norm_results.b;     %phi should be called beta for consistency
PHI_ribo = repmat(ribo_results.phi,1,6)+ribo_results.b;

%Plot NLME results and data with RNG default
tplot=0:0.5:30;
figure; grid on; hold on;
for I = 1:6
    yplot=modelfun_norm(PHI_norm(:,I)',tplot);
    plot(tplot,yplot,tplot,data,'k*','MarkerSize',1.5,'Linewidth',2)
end
xlabel("Time (h)")
ylabel("Fluorescence (a.u.)")
title("NLME results NORM model")
hold off

figure; grid on; hold on;
for I = 1:6
    yplot=modelfun_ribo(PHI_ribo(:,I)',tplot);
    plot(tplot,yplot,tplot,data,'k*','MarkerSize',1.5,'Linewidth',2) 
end
xlabel("Time (h)")
ylabel("Fluorescence (a.u.)")
title("NLME results RIBO model")
hold off



