function xdot = transfection_norm_ode(~,x,p)
xdot = zeros(2,1);

mRNA = x(1); 
GFP = x(2);

delta1 = p(1); 
pbeta = p(2); 
k2_m0_scale = p(3);

xdot(1) = - delta1*mRNA;% + dirac(t-t0);
xdot(2) = + k2_m0_scale*mRNA - pbeta*GFP;
