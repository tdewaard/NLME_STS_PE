function [texp,fluor] = extract_data(date,state,logconv)
%EXTRACT_DATA Summary of this function goes here
%   Detailed explanation goes here
switch state
    case "stable"
        state="eGFP";
    case "unstable"
        state="d2eGFP";
    otherwise
        warning("Unexpected data!")
end

disp(strcat("Loading in: ","data/",date,"_mean_",state,".xlsx"))
data=xlsread(strcat("data/",date,"_mean_",state,".xlsx"));

texp = data(:,1)/3600;
if logconv
    fluor = log(data(:,2:end));
else
    fluor = data(:,2:end);
end
end

