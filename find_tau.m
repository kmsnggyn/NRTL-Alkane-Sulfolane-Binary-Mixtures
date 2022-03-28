clear all
close all 
% clc 

set(groot,'defaultLineMarkerSize', 10, ...
    'defaultLineLineWidth', 2, ...
    'defaultAxesFontName', 'Times', ...
    'defaultAxesFontSize', 12, ...
    'defaultLegendFontSize', 12);

alkane1 = 'pentane'; lb = 0; ub = 6;
alkane2 = 'hexane';  lb = 1; ub = 9;
alkanes = {alkane1, alkane2};

for i_alkane = 1:2

alkane = alkanes{i_alkane};
fprintf(['\n\n << FOR ' upper(alkane) ' >>\n'])
load([alkane '_expdata.mat'])
load([alkane '_parameters'])

a12 = parameters_Aspen(1);
b12 = parameters_Aspen(2);
c12 = parameters_Aspen(3);
a21 = parameters_Aspen(5);
b21 = parameters_Aspen(6);
c21 = parameters_Aspen(7);

tau12_Aspen = a12 + b12./data(:,1) + c12.*log(data(:,1));
tau21_Aspen = a21 + b21./data(:,1) + c21.*log(data(:,1));

expT = data(:,1);
expX11 = data(:,2);
expX12 = data(:,3);
clear data

m = 1000;    % Adjust precision of finding optimal tau          
range = 10;

% fprintf(' m = %d\n',m)
fprintf('==========================================\n')
fprintf('  T/K      tau12     tau21        OF\n')
fprintf('------------------------------------------\n')

tic
for k = 1 : size(expT)
    
    T = expT(k);
    x11 = expX11(k);
    x12 = expX12(k);

    [tau12_opt(k), tau21_opt(k), minOF(k), OF{k}] = opt_tau(5, 2, x11, x12, m, range);
    [tau12_opt(k), tau21_opt(k), minOF(k), OF{k}] = opt_tau(tau12_opt(k), tau21_opt(k), x11, x12, m, range/m);
    [tau12_opt(k), tau21_opt(k), minOF(k), OF{k}] = opt_tau(tau12_opt(k), tau21_opt(k), x11, x12, m, range/m^2);
    [tau12_opt(k), tau21_opt(k), minOF(k), OF{k}] = opt_tau(tau12_opt(k), tau21_opt(k), x11, x12, m, range/m^3);

    fprintf(' %.2f  %f  %f  %.10f     %7f  %7f\n', ...
        T, tau12_opt(k), tau21_opt(k), minOF(k), tau12_opt(k)-tau12_Aspen(k), tau21_opt(k)-tau21_Aspen(k))
end
time = toc;
fprintf('------------------------------------------\n')
fprintf(' Elapsed time: %.2f sec\n', time)


tau12_opt = tau12_opt';
tau21_opt = tau21_opt';

fig = figure('Position',[0 10000 600 400]);
plot(expT,tau12_Aspen,'DisplayName','$\tau_{12}$, ${\mathrm {Aspen}}$'); hold on
plot(expT,tau21_Aspen,'DisplayName','$\tau_{21}$, ${\mathrm {Aspen}}$');
plot(expT,tau12_opt,'--','DisplayName','$\tau_{12}$, ${\mathrm {optimum}}$');
plot(expT,tau21_opt,'--','DisplayName','$\tau_{21}$, ${\mathrm {optimum}}$');
xlabel('$$T$$ $${\mathrm {[K]}}$$','Interpreter','latex')
ylabel('$$\tau$$', 'Interpreter','latex')
% xlim([300 390])
ylim([lb ub])
legend('Interpreter','latex')

exportgraphics(gca,[alkane '_tau_comparison.jpg'],'Resolution',300)

T_data = expT;
tau12_data = tau12_opt;
tau21_data = tau21_opt;
clear expT tau12_opt tau21_opt

opts = optimset('Display','off');
fun = @(k,xdata) k(1) + k(2)./xdata + k(3).*log(xdata) +k(4)*xdata ;

k12_0 = [6,1,0,-1];
k21_0 = [2,1,0,-1];

k12 = lsqcurvefit(fun, k12_0, T_data, tau12_data,[],[],opts);
k21 = lsqcurvefit(fun, k21_0, T_data, tau21_data,[],[],opts);
fprintf('\n')
fprintf(' Fitted coefficients:\n')
fprintf('-------------------------------------\n')
fprintf(' a12 %11.4f    ', k12(1))
fprintf(' a21 %11.4f\n', k21(1))
fprintf(' b12 %11.4f    ', k12(2))
fprintf(' b21 %11.4f\n', k21(2))
fprintf(' c12 %11.4f    ', k12(3))
fprintf(' c21 %11.4f\n', k21(3))
fprintf(' d12 %11.4f    ', k12(4))
fprintf(' d21 %11.4f\n', k21(4))
fprintf('-------------------------------------\n')

parameters_opt = [k12, k21];
save([alkane '_parameters_opt.mat'], 'parameters_opt')

tau12_fitted = fun(k12,T_data);
tau21_fitted = fun(k21,T_data);

fig = figure('Position',[0 0 600 400]);
plot(T_data,tau12_data,'x'); hold on
plot(T_data,tau21_data,'x');
plot(T_data,tau12_fitted,'-');
plot(T_data,tau21_fitted,'-');
legend(...
    '$$\tau_{12}$$, $\mathrm{optimum}$', ...
    '$$\tau_{21}$$, $\mathrm{optimum}$', ...
    '$$\tau_{12}$$, $\mathrm{fitted}$', ...
    '$$\tau_{21}$$, $\mathrm{fitted}$', ...
    'Interpreter', 'latex','Location','northeast')
xlabel('$$T$$ $${\mathrm {[K]}}$$','Interpreter','latex')
ylabel('$$\tau$$', 'Interpreter','latex')
% xlim([300 390])
ylim([lb ub])
exportgraphics(gca,[alkane '_tau_fitted.jpg'],'Resolution',300)

end



function [tau12_opt, tau21_opt, minOF2, OF] = opt_tau(tau12_in, tau21_in, x1alpha, x1beta, m, range)

    alpha = 0.3;

    for i = 1: m
        tau12 = tau12_in + range * (-0.5 + i/m);

        for j = 1: m
            tau21 = tau21_in + range * (-0.5 + j/m);

            G12 = exp(-alpha*tau12);
            G21 = exp(-alpha*tau21);
        
            x1 = x1alpha;
            x2 = 1 - x1;
            gamma12 = exp( x2^2 * ( tau21*( G21/(x1+x2*G21) )^2 + G12*tau12 / (x2+x1*G12)^2 ) );
            gamma22 = exp( x1^2 * ( tau12*( G12/(x2+x1*G12) )^2 + G21*tau21 / (x1+x2*G21)^2 ) );
            
            k11 = x1*gamma12;
            k21 = x2*gamma22;
            
            x1 = x1beta;
            x2 = 1 - x1;
            gamma11 = exp( x2^2 * ( tau21*( G21/(x1+x2*G21) )^2 + G12*tau12 / (x2+x1*G12)^2 ) );
            gamma21 = exp( x1^2 * ( tau12*( G12/(x2+x1*G12) )^2 + G21*tau21 / (x1+x2*G21)^2 ) );
            
            k12 = x1*gamma11;
            k22 = x2*gamma21;
        
            OF(i,j) = (log(k11) - log(k12))^2 + (log(k21) - log(k22))^2;
            
        end
    end

    [minOF1, index1] = min(OF);
    [minOF2, index2] = min(minOF1);
    ind1 = index1(index2);
    ind2 = index2;
    tau12_opt = tau12_in + range * (-0.5 + ind1/m);
    tau21_opt = tau21_in + range * (-0.5 + ind2/m);
end