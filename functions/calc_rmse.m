function [rmse] = calc_rmse(data,texp,y)
%CALC_RMSE Summary of this function goes here
%   Calculates RMSE for given model output y and data cobination
rmse=sqrt(sum((y-data).^2)/length(texp));
end

