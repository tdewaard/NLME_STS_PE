%% Loading in data & initialise models
clear; close all
addpath('functions')

dates = ["20160427","20160513", "20160613"];

[texp,fluor] = extract_data(dates(2),'stable',true);

MODELS = init_models();
simulate_model(MODELS(2),texp,[],MODELS(2).p0,true,'r',3);

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
