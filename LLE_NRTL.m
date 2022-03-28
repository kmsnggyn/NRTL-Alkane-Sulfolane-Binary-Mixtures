%% Full LLE Diagram

clear all
close all 
clc 

alkanes = {'pentane', 'hexane'};
Tlims = [300, 400; 300, 480];

for i_alkane = 2:2

alkane = alkanes{i_alkane};
Tlb = Tlims(i_alkane,1);  Tub = Tlims(i_alkane,2);

load([alkane '_expdata.mat'])
load([alkane '_tau_Aspen.mat'])
load([alkane '_parameters'])
load([alkane '_parameters_opt.mat'])

set(groot,'defaultLineMarkerSize', 10, ...
    'defaultLineLineWidth', 2, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesFontSize', 12);

T_data = data(:,1);
x11_data = data(:,2);
x12_data = data(:,3);


tic

for par = 1:2
clear result

    fprintf('\n\n')
    if par ==1
        parameters = parameters_Aspen;
        fprintf(' Parameters from Aspen')
    elseif par ==2
        parameters = parameters_opt;
        fprintf(' Parameters optimized')
    end

fprintf('\n')
fprintf('=======================\n')
fprintf('  T/K     x11     x12\n')
fprintf('-----------------------\n')

x12_in = 0;
x11_in = 1;

for k = 1:2000
    
    T = 300 + 0.1*(k-1) ;
    result(k,1) = T;
    [result(k,3), result(k,2), ~, ~, ~] = NRTL(x12_in, x11_in, parameters, T);

    if result(k,3) >= result(k,2)
        fprintf('-----------------------\n')
        result(k,:) = [];
        k = k-1;
        break
    end
    
    if ((k-1)/ 20) == fix((k-1) / 20)
        fprintf(' %.1f  %0.4f  %0.4f\n', T, result(k,2), result(k,3))
    end
    x12_in = result(k,3);
    x11_in = result(k,2);
end
fprintf('\n')
toc

LLE = figure('Position',[0 0 500 500]);

plot(x12_data, T_data, 'x'); hold on
plot(x11_data, T_data, 'x');
plot(result(:,3), result(:,1), '-');
plot(result(:,2), result(:,1), '-');
pbaspect([1 1 1])
xlim([0 1])
ylim([Tlb Tub])
legend(...
    '$$x_1^\alpha$$, $${\mathrm {exp}}$$', ...
    '$$x_1^\beta$$, $${\mathrm {exp}}$$', ...
    '$$x_1^\alpha$$, $${\mathrm {calc}}$$', ...
    '$$x_1^\beta$$, $${\mathrm {calc}}$$ ', ...
    'Interpreter','latex', ...
    'Location','northwest')
xlabel('$$x_1$$','Interpreter','latex')
ylabel('$$T$$ $${\mathrm {[K]}}$$','Interpreter','latex')
exportgraphics(gca,[alkane '_LLE_' num2str(par) '.jpg'],'Resolution',300)
end
end
%% NRTL Result for Error Analysis

clear all
close all 
% clc 

% alkane = 'pentane'; Tlb = 300;  Tub = 400;
alkane = 'hexane';  Tlb = 300;  Tub = 480;

load([alkane '_expdata.mat'])
load([alkane '_tau_Aspen.mat'])
load([alkane '_parameters'])
load([alkane '_parameters_opt.mat'])

set(groot,'defaultLineMarkerSize', 10, ...
    'defaultLineLineWidth', 2, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesFontSize', 12);

T_data = data(:,1);
x11_data = data(:,2);
x12_data = data(:,3);

for par = 1:2
    fprintf('\n\n')
    if par ==1
        parameters = parameters_Aspen;
        fprintf([' < ' upper(alkane) '/SULFOLANE > w/ parameters from literature'])
    else
        parameters = parameters_opt;
        fprintf([' < ' upper(alkane) '/SULFOLANE > w/ parameters optimized'])
    end
    fprintf('\n')
    fprintf('============================================================\n')
    fprintf('           Exp     NRTL              Exp     NRTL\n')
    fprintf(' T/K       x11     x11     err(%%)    x12     x12     err(%%)\n')
    fprintf('------------------------------------------------------------\n')
    
    x12_in = 0;
    x11_in = 1;
    tic
    
    for k = 1:size(T_data)
            
        T = T_data(k);
        result(k,1) = T;
    
        [result(k,3), result(k,2), tau12(k), tau21(k), OF(k)] = NRTL(x12_in, x11_in, parameters, T);
        err(k,1) = abs( ( result(k,2) - data(k,2) ) / data(k,2) ) *100 ; % x11 error
        err(k,2) = abs( ( result(k,3) - data(k,3) ) / data(k,3) ) *100 ; % x12 error
        x12_in = result(k,3);
        x11_in = result(k,2);
        fprintf(' %.2f    %.4f  %.4f  %.4f    %.4f  %.4f  %0.4f\n', ...
            T, x11_data(k), result(k,2), err(k,1), x12_data(k), result(k,3), err(k,2))
    
    end
    fprintf('------------------------------------------------------------\n')

    time = toc;
    meanerr(1) = mean(err(:,1));
    meanerr(2) = mean(err(:,2));
    toterr = mean(meanerr);
    
    fprintf('     x11 error %7.4f\n', meanerr(1))
    fprintf('     x12 error %7.4f\n', meanerr(2))
    fprintf(' overall error %7.4f\n', toterr)
    fprintf('  elapsed time %7.4f seconds\n', time)
    


end

%% Function

function [x12_cal, x11_cal, tau12, tau21, minOF2] = NRTL(x12_in,x11_in,parameters, T)

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

n = 10*10^4;
x12 = zeros(n,1);
x11 = zeros(n,1);
% OF = zeros(n/2);
stop_i = 0;
stop_j = 0;

for i = 1:n
    x12(i) = x12_in + (i-1)/n;

    for j = 1:n
        if stop_j ==1
            stop_j = 0;
            break
        end
        x11(j) = x11_in - (j-1)/n;

        % alpha phase (rich in 2)
        x1 = x12(i);
        x2 = 1 - x1;
        gamma11 = exp( x2^2 * ( tau21*( G21/(x1+x2*G21) )^2 + G12*tau12 / (x2+x1*G12)^2 ) );
        gamma21 = exp( x1^2 * ( tau12*( G12/(x2+x1*G12) )^2 + G21*tau21 / (x1+x2*G21)^2 ) );
        
        k11 = x1*gamma11;
        k21 = x2*gamma21;
        
        % beta phase (rich in 1)
        x1 = x11(j);
        x2 = 1 - x1;
        gamma12 = exp( x2^2 * ( tau21*( G21/(x1+x2*G21) )^2 + G12*tau12 / (x2+x1*G12)^2 ) );
        gamma22 = exp( x1^2 * ( tau12*( G12/(x2+x1*G12) )^2 + G21*tau21 / (x1+x2*G21)^2 ) );
        
        k12 = x1*gamma12;
        k22 = x2*gamma22;

        OF(i,j) = (k11 - k12)^2 + (k21 - k22)^2;

        if j>2
            if OF(i,j) > OF(i,j-1)
                minOF(i) = OF(i,j-1);
                min_j = j-1;
                stop_j = 1;
            end
        end
    end
%     fprintf('\n%f', min(OF(i)));

    if i> 2
        if minOF(i) > minOF(i-1)
            min = minOF(i-1);
            min_i = i-1;
            stop_i = 1;
        end
    end

    if stop_i == 1
        break
    end
end

minOF2 =  min;
x11_cal = x11(min_j);
x12_cal = x12(min_i);            

end
