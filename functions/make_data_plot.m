function [h1,h2] = make_data_plot(data,time,conf)
%MAKE_DATA_PLOT Summary of this function goes here
%   Detailed explanation goes here
if (~exist('conf', 'var'))
        conf = false;
end

figure;

for k = 1:length(data) %for-loop to plot with transparent lines, set in lh.Color
    lh=plot(time,data(:,k),'k'); hold on
    lh.Color=[0,0,0,0.1];   %black, 4th number sets transparency
end

set(findall(gcf,'-property','FontSize'),'FontSize',13)
title("\bf{Experimental Fluorescence Time Profiles}")
xlim([0,time(end)])
ylim([0,max(data(:))+1])
xlabel('Time (h)')
ylabel('Fluorescence intensity (a.u.)')

if conf
    disp('Plotting mean and standard deviation...')
    h1=plot(time,mean(data'),'r', 'LineWidth',2);
    h2=plot(time,mean(data')+std(data'),'r--', time, mean(data')-std(data'),'r--', 'LineWidth',1);
    legend([h1,h2(1)],'average','stdev')
end

hold off
end

