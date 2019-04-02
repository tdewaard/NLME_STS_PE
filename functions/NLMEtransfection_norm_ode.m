function xdot = NLMEtransfection_norm_ode(~,x,p)
xdot = zeros(2,1);

%p = [delta1,pbeta,k2_m0_scale,t0,offset];

mRNA = x(1); 
GFP = x(2);

delta1 = p(:,1); 
pbeta = p(:,2); 
k2_m0_scale = p(:,3);

xdot(1) = - delta1*mRNA;% + dirac(t-t0);
xdot(2) = + k2_m0_scale*mRNA - pbeta*GFP;
