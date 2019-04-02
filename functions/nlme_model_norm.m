function [y] = nlme_model_norm(PHI,texp,logconv)
%NLME_MODEL Summary of this function goes here
%   Detailed explanation goes here

%p = [delta1,pbeta,k2_m0_scale,t0,offset];

if logconv
    t0 = exp(PHI(:,4));
    x0 = exp([1; 0]);
    offset= exp(PHI(:,5));

    [t_out,x] = ode15s(@NLMEtransfection_norm_ode,[t0,texp(end)], x0, [], exp(PHI));
else
    t0 = PHI(:,4);
    x0 = [1; 0];
    offset= PHI(:,5);

    [t_out,x] = ode15s(@NLMEtransfection_norm_ode,[t0,texp(end)], x0, [], PHI);
end

deltat=texp(2)-texp(1);
t_out=[0; t_out(1)-deltat; t_out];  %make time vector start at 0
x2 = [0; 0; (x(:,2))];    %zero padding of output for the time until t0

%output

y = log(x2 + offset);
y = interp1(t_out,y,texp); %interpolate such that simulated y and data match

if any(isnan(y))
    if any(~isfinite(y))
        disp("Y has Inf and NaN values!")
    else
        disp("Y has NaN values!")
    end
end
%     y = real(y);
%     y(isnan(y)|~isfinite(y))=0;

end
