%% Loading in data & initialise models
clear; close all
addpath('functions')

dates = ["20160427","20160513", "20160613"];

[texp,fluor] = extract_data(dates(2),'stable',true);

%convert real data to 'grouped data' format for NLMEfit
TIME=repmat(texp,length(fluor),1)';
TIME=TIME(:);
GROUPS=repmat(1:length(fluor),length(texp),1);
GROUPS=GROUPS(:);
DATA=fluor(:);

MODELS = init_models();     %MODELS(1) = NORM, MODELS(2) = RIBO MODELS(3) = ENZDEGRIBO.
simulate_model(MODELS(3),texp,[],MODELS(3).p0,true,'r',3);

%% NLME on real data

n_groups=2;
index=length(texp)*n_groups;
t=TIME(1:index);
d=DATA(1:index);
g=GROUPS(1:index);
m=MODELS(2);  
filename=strcat("Results\RealData",num2str(n_groups),"G_Result_NLMEFIT_",m.name,".mat")
[Beta,Psi,stats,b,t_elapsed]=nlme_param_est(m,t,d,g,false)
save(filename,'Beta','Psi','stats','b','t_elapsed')

%% NLME on synth data, different group sizes

%for RIBO it only works until n_groups = 5???
m=MODELS(2);
n_groups=3;
sigma_xi=0.05;

[sTIME, sGROUPS, sDATA]=generate_synth_data(m,sigma_xi,1,false,n_groups,true);
%%
%NLMEfit on synthetic data
[Beta,Psi,stats,b,t_elapsed]=nlme_param_est(m,sTIME(:),sDATA(:),sGROUPS(:),false)
filename=strcat("SynthData",num2str(n_groups),"G_xi=",num2str(sigma_xi),"Result_NLMEFIT_",m.name,".mat");
save(filename,'Beta','Psi','stats','b','t_elapsed')

