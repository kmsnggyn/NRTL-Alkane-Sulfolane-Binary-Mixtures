clear all
close all
clc 

set(groot,'defaultLineMarkerSize', 10, ...
    'defaultLineLineWidth', 1, ...
    'defaultAxesFontName', 'Times', ...
    'defaultAxesFontSize', 12);

alkane = 'pentane'; T_range = [300 325 350 375 400 425];
% alkane = 'hexane';  T_range = [380 400 420 440 460 480];

load([alkane '_expdata.mat'])
load([alkane '_parameters'])
load([alkane '_parameters_opt.mat'])

parameters = parameters_Aspen;
parameters = parameters_opt;


expT = data(:,1);
expX11 = data(:,2);
expX12 = data(:,3);

Gxplot = figure('Position',[0 0 400 600]);
dGxplot = figure('Position',[400 0 400 600]);
ddGxplot = figure('Position',[800 0 400 600]);

figure(2),plot(linspace(0,1,1001),0,'-k');
m = 1000000;
stop = 0;
for j = 1:6
    if stop == 1
        break
    end
    T = T_range(j);

    for i = 1:m+1
        x1(i) = (i-1)/m;
        deltaGRT(i) = G_NRTL(parameters, T,x1(i));
    end

    for i = 1:m
        d_deltaGRT(i) = ( deltaGRT(i+1) - deltaGRT(i) ) / (1/m);
    end

    for i = 1:m-1
        dd_deltaGRT(i) = ( d_deltaGRT(i+1) - d_deltaGRT(i) ) / (1/m);
    end
    
    figure(1),plot(x1,deltaGRT,'DisplayName',[num2str(T) 'K']); hold on
    xlim([0 1]);
    figure(2),plot(x1(1:m),d_deltaGRT,'DisplayName',[num2str(T) 'K']); hold on
    xlim([0 1]); ylim([-5 5])
    figure(3),plot(x1(1:m-1),dd_deltaGRT,'DisplayName',[num2str(T) 'K']); hold on
    xlim([0 1]); ylim([-5 5])
    
end

figure(1),legend('Location','northwest','AutoUpdate','off')
xlabel('$$x_1$$','Interpreter','latex')
exportgraphics(gca,[alkane '_G.jpg'],'Resolution',300)

figure(2),legend('Location','southwest','AutoUpdate','off')
yline(0,'LineWidth',1)
xlabel('$$x_1$$','Interpreter','latex')
exportgraphics(gca,[alkane '_dG.jpg'],'Resolution',300)

figure(3),legend('Location','southwest','AutoUpdate','off')
yline(0,'LineWidth',1)
xlabel('$$x_1$$','Interpreter','latex')
exportgraphics(gca,[alkane '_ddG.jpg'],'Resolution',300)

function delta_GRT = G_NRTL(parameters, T,x1)

alpha = 0.3;
a12 = parameters(1);
b12 = parameters(2);
c12 = parameters(3);
d12 = parameters(4);
a21 = parameters(5);
b21 = parameters(6);
c21 = parameters(7);
d21 = parameters(8);

tau12 = a12 + b12/T + c12*log(T) + d12*T;
tau21 = a21 + b21/T + c21*log(T) + d21*T;
G12 = exp(-alpha*tau12);
G21 = exp(-alpha*tau21);

x2 = 1 - x1;
delta_GRT = ( x1*log(x1) + x2*log(x2) ...
    + x1*x2*( (G21*tau21)/(x1+x2*G21) + (G21*tau12)/(x2+x1*G12) )); 

end