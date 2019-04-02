function [MODELS] = init_models()
%INIT_MODELS Summary of this function goes here
%   Detailed explanation goes here
%Params
t0 = 1; 
delta1=1e-1;
pbeta=1e-1;
k2_m0_scale=1e3;
offset=10;
k1_m0=1;
frac_R0_m0=1;
k2=1;

delta1_m0=1e-3;
frac_E0_m0=0.5;
delta3=1e-3;

%Initialise model objects
pnames={'\delta','k_{1}m_{0}','k_2','^{r_{0}}/_{m_{0}}','k_{2}m_{0}scale', '\gamma', 'offset', 't_0'};

NORM=Model;
NORM.name='NORM';
NORM.pnames={pnames{1},pnames{6},pnames{5},pnames{8},pnames{7}};
NORM.p0=[delta1,pbeta,k2_m0_scale,t0,offset];
NORM.odehandle=@transfection_norm_ode;
NORM.errorhandle=@transfection_norm_error;
NORM.nlmehandle=@(PHI,t)nlme_model_norm(PHI,t,true);

RIBO=Model;
RIBO.name='RIBO';
RIBO.pnames=pnames;
RIBO.p0=[delta1, k1_m0, k2, frac_R0_m0, k2_m0_scale, pbeta, offset, t0];
RIBO.odehandle=@transfection_ribo_ode;
RIBO.errorhandle=@transfection_ribo_error;
RIBO.nlmehandle=@(PHI,t)nlme_model_ribo(PHI,t,true,false);

ENZDEGRIBO=Model;
ENZDEGRIBO.name='ENZDEGRIBO';
ENZDEGRIBO.pnames={'\delta_{1}m_{0}','\gamma','k_{2}m_{0}scale','t_0','offset','k_{1}m_{0}','^{r_{0}}/_{m_{0}}','k_2','^{E_{0}}/_{m_{0}}','\delta_3'};
ENZDEGRIBO.p0=[delta1_m0,pbeta,k2_m0_scale,t0,offset,k1_m0,frac_R0_m0,k2,frac_E0_m0,delta3];
ENZDEGRIBO.odehandle=@transfection_enzdegribo_ode;
ENZDEGRIBO.errorhandle=[];
ENZDEGRIBO.nlmehandle=@(PHI,t)nlme_model_enzdegribo(PHI,t,true,false);

MODELS=[NORM RIBO ENZDEGRIBO];
end

