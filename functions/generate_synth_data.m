function [sTIME,sGROUPS, sDATA] = generate_synth_data(model,sigma_xi,p_var,random,ngroups,makeplot)
%GENERATE_SYNTH_DATA Summary of this function goes here
%   Generates synthetic data with RIBO model. sTIME, sGROUPS and sDATA are
%   column vectors ready for implementation in nlmefit()
%   sigma_xi = standard deviation of added noise
%   p_var = coefficient that adds more variation between individuals' model
%   parameters
if ~random    
    rng default
end

phi0=zeros(ngroups,length(model.p0));
switch model.name
    case 'NORM'
        %p0 = [delta1,pbeta,k2_m0_scale,t0,offset];
        for i=1:ngroups
            phi0(i,:)=[normrnd(model.p0(1),p_var*1e-2),normrnd(model.p0(2),p_var*0.01),normrnd(model.p0(3),p_var*100),normrnd(model.p0(4),p_var*0.01), normrnd(model.p0(5),p_var*2)];
        end
    case 'RIBO'
        %p0=[delta1, k1_m0, k2, frac_R0_m0, k2_m0_scale, pbeta, offset, t0];
        for i=1:ngroups
            phi0(i,:)=[normrnd(model.p0(1),p_var*1e-2), normrnd(model.p0(2),p_var*0.05), normrnd(model.p0(3),p_var*0.2), normrnd(model.p0(4),p_var*0.5),...
                       normrnd(model.p0(5),p_var*100), normrnd(model.p0(6),p_var*0.01), normrnd(model.p0(7),p_var*2), normrnd(model.p0(8),p_var*0.01)];
        end
end


t=0:0.5:30;
y=zeros(ngroups,length(t));
data=zeros(ngroups,length(t));
for i=1:ngroups
    y(i,:)=model.nlmehandle(log(phi0(i,:)),t);
    xi=normrnd(0,sigma_xi,[1,length(t)]);
    data(i,:)=y(i,:)+xi;
end

TIME=repmat(t,ngroups,1)';
sTIME=TIME;
GROUPS=repmat(1:ngroups,length(t),1);
sGROUPS=GROUPS;
DATA=data';
sDATA=DATA;
disp(strcat("Data generate succesfully with ",num2str(ngroups)," groups and sigma_xi = ",num2str(sigma_xi),"!"))

if makeplot
    plot(t,data,'Color',[0 0 0 0.5])
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    title(strcat("\bf{Synthetic data ",num2str(ngroups)," groups, \sigma_{\xi} = ",num2str(sigma_xi),"}"))
    xlabel("Time (h)")
    ylabel("Fluorescence (a.u.)")
    hold off
end


end

