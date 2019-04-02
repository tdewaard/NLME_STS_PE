%transfection_norm_error
function e = transfection_norm_error(p, texp, data)

%p = [delta1,pbeta,k2_m0_scale,t0,offset];
t0 = p(4);
offset = p(5);

%x = [mRNA,GFP];
x0 = [1 0];


[t,x] = ode15s(@transfection_norm_ode,[t0 texp(end)], x0, [], p);

deltat=texp(2)-texp(1);
t=[0; t(1)-deltat; t];  %make time vector start at 0
x2 = [0; 0; x(:,2)];    %zero padding of output for the time until t0

%output

y = log(x2 + offset);
y = interp1(t,y,texp); %interpolate such that simulated y and data match
e = y - data;
end
