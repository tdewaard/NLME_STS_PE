function err = transfection_ribo_error(p, texp, data)
%ERR Summary of this function goes here
%   Detailed explanation goes here
%
%p0=[delta1, k1_m0, k2, frac_R0_m0, k2_m0_scale, pbeta, offset, t0];
t0=p(8);
offset=p(7);
x0=[1; 0; p(4)];

[t,x] = ode15s(@transfection_ribo_ode,[t0 texp(end)], x0, [], p);
deltat=texp(2)-texp(1);
t=[0; t(1)-deltat; t];  %make time vector start at 0
x2 = [0; 0; (x(:,2))];    %zero padding of output for the time until t0

%output
y = log(x2 + offset);
y = interp1(t,y,texp); %interpolate such that simulated y and data match
err = y - data;
end

