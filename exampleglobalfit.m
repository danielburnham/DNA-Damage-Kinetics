clear all; close all; clc;
addpath(genpath('/Users/danielburnham/Dropbox/anal/altmany-export_fig-04ca93c'))

load('raw50ng.mat');
load('stdev50ng.mat');

raw = raw50ng;
stdev = stdev50ng;

G_data = raw(:,1);
S_data = raw(:,2)/100;
R_data = raw(:,3)/100;
L_data = raw(:,4)/100;

S_err = stdev(:,2)/100;
R_err = stdev(:,3)/100;
L_err = stdev(:,4)/100;

S0 = 1;
% R0 = 0;
% L0 = 0;

G = linspace(0,40,1001);
% options = optimset('MaxFunEvals',100,'MaxIter',100,'disp','iter');
options = optimset('display','iter');

%% global
[p_global_nl,R,J,covB,~,~] = nlinfit([G_data; G_data; G_data],[S_data; R_data; L_data],@(param,G_data)(slr_global_nl(param,G_data)),[1 0.001 1],options,'weights',[(1./(S_err.^2)); (1./(R_err.^2)); (1./(L_err.^2))])
p_global_nl
p_global_nl_SE = diag(sqrt(covB))
ci = nlparci(p_global_nl,R,'jacobian',J)
ci2 = nlparci(p_global_nl,R,'covar',covB)

S_nl = p_global_nl(3)*exp(-p_global_nl(1)*G);
R_nl = (p_global_nl(1)/(p_global_nl(1)-p_global_nl(2)))*(exp(-p_global_nl(2)*G)-exp(-p_global_nl(1)*G))*p_global_nl(3);
% R_nl_lower = (ci(1,2)/(ci(1,2)-ci(2,2)))*(exp(-ci(2,2)*G)-exp(-ci(1,2)*G))*S0;
L_nl = (1+(((p_global_nl(1)*exp(-p_global_nl(2)*G))-(p_global_nl(2)*exp(-p_global_nl(1)*G)))/(p_global_nl(2)-p_global_nl(1))))*p_global_nl(3);
% % (1+(((param2(2)*exp(-param2(2)*G_data))-(param2(2)*exp(-param2(2)*G_data)))/(param2(1)-param2(2))))*S0
% [p_global_nl2,R2,J2,covB2,~,~] = nlinfit(G_data,L_data,@(param2,G_data)( (1+(((param2(1)*exp(-param2(2)*G_data))-(param2(2)*exp(-param2(1)*G_data)))/(param2(2)-param2(1))))*S0 ),[1.5,0.02],'weights',[(1./(L_err.^2))]);
% p_global_nl2
% p_global_nl_SE2 = diag(sqrt(covB2))
% ci2 = nlparci(p_global_nl2,R2,'jacobian',J2)
% ci22 = nlparci(p_global_nl2,R2,'covar',covB2)

%% plots
%% data
f1 = figure(1);
errorbar(G_data,S_data,S_err,'bo','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')
hold on
errorbar(G_data,R_data,R_err,'rs','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')
errorbar(G_data,L_data,L_err,'k^','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')

set(gca,'linewidth',1.4,'FontName','arial','fontsize',12)
axis([-2 42 -0.04 1.04])
xlabel('Dose (Gy)','FontName','arial','fontsize',14);
ylabel('Conc (%/100)','FontName','arial','fontsize',14);

leg = legend('supercoiled','relaxed','linear','Location','northeast','Orientation','vertical');
leg.FontName = 'arial';
leg.FontSize = 14;
legend('boxoff')
box on
set(gca, 'Layer', 'top')

set(gcf,'Color','w')
export_fig srl_data_50ng.pdf -painters -rgb


%% linear fit
f2 = figure(2);
errorbar(G_data,S_data,S_err,'bo','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')
hold on
errorbar(G_data,R_data,R_err,'rs','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')
errorbar(G_data,L_data,L_err,'k^','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')

plot(G,S_nl,'b-','linewidth',1.2)
plot(G,R_nl,'r-','linewidth',1.2)
plot(G,L_nl,'k-','linewidth',1.2)

% set(gca,'Position',[0.13 0.11 0.55 0.815]);

set(gca,'linewidth',1.4,'FontName','arial','fontsize',12)
axis([-2 42 -0.04 1.04])
xlabel('Dose (Gy)','FontName','arial','fontsize',14);
ylabel('Conc (%/100)','FontName','arial','fontsize',14);

leg = legend('supercoiled data','relaxed data','linear data','supercoiled fit','relaxed fit','linear fit','Location','northeast','Orientation','vertical');
leg.FontName = 'arial';
leg.FontSize = 14;
legend('boxoff')
box on
set(gca, 'Layer', 'top')

to_add_1 = strcat(['k_{sr} = ' num2str(p_global_nl(1),'%0.2f') '\pm' num2str(p_global_nl_SE(1),'%0.1g') ' Gy^{-1}']);
text(25,0.3,to_add_1,'FontName','arial','fontsize',14);
to_add_2 = strcat(['k_{rl} = ' num2str(p_global_nl(2),'%0.3f') '\pm' num2str(p_global_nl_SE(2),'%0.1g') ' Gy^{-1}']);
text(25,0.2,to_add_2,'FontName','arial','fontsize',14);

set(gcf,'Color','w')
export_fig srl_global_fit_linear_50ng.pdf -painters -rgb

%% log fit
f3 = figure(3);
errorbar(G_data,S_data,S_err,'bo','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')
hold on
errorbar(G_data,R_data,R_err,'rs','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')
errorbar(G_data,L_data,L_err,'k^','linewidth',1.2,'MarkerSize',12,'MarkerFaceColor','w')

plot(G,S_nl,'b-','linewidth',1.2)
plot(G,R_nl,'r-','linewidth',1.2)
% plot(G,R_nl_lower,'r--','linewidth',1.2)
plot(G,L_nl,'k-','linewidth',1.2)

set(gca,'Position',[0.13 0.11 0.55 0.815]);

set(gca,'linewidth',1.4,'FontName','arial','fontsize',12)
axis([0.05 50 -0.04 1.04])
xlabel('Dose (Gy)','FontName','arial','fontsize',14);
ylabel('Conc (%/100)','FontName','arial','fontsize',14);

leg = legend('supercoiled data','relaxed data','linear data','supercoiled fit','relaxed fit','linear fit','Location','northeast','Orientation','vertical');
leg.FontName = 'arial';
leg.FontSize = 14;
legend('boxoff')
rect = [0.70, 0.65, .25, .25];
set(leg, 'Position', rect)
box on
set(gca, 'Layer', 'top','XScale','log')

to_add_1 = strcat(['k_{sr} = ' num2str(p_global_nl(1),'%0.2f') '\pm' num2str(p_global_nl_SE(1),'%0.1g') ' Gy^{-1}']);
text(90,0.3,to_add_1,'FontName','arial','fontsize',14);
to_add_2 = strcat(['k_{rl} = ' num2str(p_global_nl(2),'%0.3f') '\pm' num2str(p_global_nl_SE(2),'%0.1g') ' Gy^{-1}']);
text(90,0.2,to_add_2,'FontName','arial','fontsize',14);

set(gcf,'Color','w')
export_fig srl_global_fit_log_50ng.pdf -painters -rgb

function out_nl = slr_global_nl(x,dose)
part = length(dose)/3;
S = x(3)*exp(-x(1)*dose(1:part));
R = (x(1)/(x(2)-x(1)))*(exp(-x(1)*dose(part+1:(2*part)))-exp(-x(2)*dose(part+1:(2*part))))*x(3);
L = (1+(((x(1)*exp(-x(2)*dose((2*part)+1:end)))-(x(2)*exp(-x(1)*dose((2*part)+1:end))))/(x(2 )-x(1))))*x(3);
out_nl = [S; R; L];
end