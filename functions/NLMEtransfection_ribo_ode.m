function [xdot] = NLMEtransfection_ribo_ode(~,y, p)
%Fröhlich model ii) TRANSFORMED ODEs

%p0=[delta1, k1_m0, k2, frac_R0_m0, k2_m0_scale, pbeta, offset, t0];

delta1=p(:,1);
k1_m0=p(:,2);
k2=p(:,3);
frac_R0_m0=p(:,4);
k2_m0_scale=p(:,5);
pbeta=p(:,6);

MRNA=y(1);
GFP=y(2);
RIBO=y(3);

xdot=zeros(3,1);
xdot(1) = - delta1*MRNA - k1_m0*MRNA*RIBO + k2*(frac_R0_m0 - RIBO);                 %mrna ode
xdot(2) = + k2_m0_scale*(frac_R0_m0 - RIBO) - pbeta*GFP;                            %gfp ode, pbeta=gamma?
xdot(3) = - k1_m0*MRNA*RIBO + k2*(frac_R0_m0 - RIBO);                               %ribo ode
end

