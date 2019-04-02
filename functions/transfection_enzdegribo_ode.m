function [xdot] = transfection_enzdegribo_ode(~,x,p)
%TRANSFECTION_ENZDEGRIBO_ODE Summary of this function goes here
%   Detailed explanation goes here

%x = [MRNA,GFP,ENZ,RIBO];
MRNA = x(1); 
GFP = x(2);
ENZ = x(3);
RIBO = x(4);

%p = [delta1_m0,pbeta,k2_m0_scale,t0,offset,k1_m0,frac_R0_m0,k2,frac_E0_m0,delta3];
delta1_m0 = p(1); 
pbeta = p(2); 
k2_m0_scale = p(3);
%t0=p(:,4);
%offset = p(:,5);
k1_m0 = p(6);
frac_R0_m0 = p(7);
k2 = p(8);
frac_E0_m0 = p(9);
delta3 = p(10);

xdot=zeros(4,1);
xdot(1) = - delta1_m0*MRNA*ENZ - k1_m0*MRNA*RIBO + k2*(frac_R0_m0-RIBO);% + dirac(t-t0);% + 1/(st*sqrt(2*pi))*exp(-(t-t0)^2/(2*st^2));
xdot(2) = + k2_m0_scale*(frac_R0_m0-RIBO) - pbeta*GFP;
xdot(3) = - delta1_m0*MRNA*ENZ + delta3*(frac_E0_m0-ENZ);
xdot(4) = - k1_m0*MRNA*RIBO + k2*(frac_R0_m0-RIBO);
end

