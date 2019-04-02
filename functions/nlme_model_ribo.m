function [y] = nlme_model_ribo(PHI,texp,logconv,nlme_error)
%NLME_MODEL Summary of this function goes here
%   Detailed explanation goes here

%p0=[delta1, k1_m0, k2, frac_R0_m0, k2_m0_scale, pbeta, offset, t0];

if logconv
    t0 = exp(PHI(:,8));
    x0 = exp([1; 0; PHI(:,4)]);
    offset= exp(PHI(:,7));

    [t_out,x] = ode15s(@NLMEtransfection_ribo_ode,[t0,texp(end)], x0, [], exp(PHI));
else
    t0 = PHI(:,8);
    x0 = [1; 0; PHI(:,4)];
    offset= PHI(:,7);

    [t_out,x] = ode15s(@NLMEtransfection_ribo_ode,[t0,texp(end)], x0, [], PHI);
end

deltat=texp(2)-texp(1);
t_out=[0; t_out(1)-deltat; t_out];  %make time vector start at 0
x2 = [0; 0; (x(:,2))];    %zero padding of output for the time until t0

%output

y = log(x2 + offset);
y = interp1(t_out,y,texp);
if any(isnan(y))
    if any(~isfinite(y))
        disp("Y has Inf and NaN values!")
    else
        disp("Y has NaN values!")
    end
end
if nlme_error
    y = real(y);
    y(isnan(y)|~isfinite(y))=0;
end

