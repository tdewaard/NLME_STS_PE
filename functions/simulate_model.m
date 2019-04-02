function [t,x,y,h] = simulate_model(model,texp,opts,params,make_plot,color,linewidth)
%SIMULATE_MODEL Summary of this function goes here
%   Detailed explanation goes here
switch model.name
    case 'NORM'
        t0 = params(4);
        x0 = [1; 0];
        offset=params(5);
    case 'RIBO'
        t0 = params(8);
        x0 = [1; 0; params(4)];
        offset=params(7); 
    case 'ENZDEGRIBO'
        t0=params(4);
        offset = params(5);
        x0=[1; 0; params(9); params(7)];
end


[t,x] = ode15s(model.odehandle,[t0, texp(end)], x0, opts, params);
deltat=texp(2)-texp(1);
t=[0; t(1)-deltat; t];  %make time vector start at 0
x2 = [0; 0; x(:,2)];    %zero padding of output for the time until t0

%output

y = log(x2 + offset);
y = interp1(t,y,texp); %interpolate such that simulated y and data match

if make_plot
    h=plot(texp,y,'LineWidth',linewidth,'color',color);
    xlabel('Time (h)')
    ylabel('Fluorescence intensity (a.u.)')
end
end

