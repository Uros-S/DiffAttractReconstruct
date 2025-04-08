% Uros Sutulovic, 03/2025

clear; close all; clc;
addpath([pwd,'/Systems/'])
addpath([pwd,'/Util/'])
figID = 0;

noise_type = 5;         % 1 = harmonic
                        % 2 = unbounded
                        % 3 = white Gaussian variance 0.01
                        % 4 = white Gaussian variance 0.1
                        % 5 = white Gaussian variance 1

additive_noise = 1;     % 0 = multiplicative noise
                        % 1 = additive noise

% Differentiator parameters
n_d = 2;
n_f = 10;
L_0 = 50;
q = 0;              % NaN to disable discrete L-adaptation
trans_diff = 5;     % differentiator transient to discard 

% Simulation time
dt = 1e-3;
t0 = 0;
tf = 3000;
T = length(0:dt:tf);                % number of data points
st = 1e-3;                           % sampling time
N = length(0:round(st/dt):tf/dt);   % sampled data points

%% Epileptor model simulation
% Initial conditions
x(:,1) = [0.022;0.91;3.84;-1.11;0.73;0];
taux(1) = t0;

% Parameters
x0 = -1.6;      y0 = 1;
Irest1 = 3.1;   Irest2 = 0.45;
tau0 = 2857;    tau2 = 10;
gamma = 0.01;

% Simulation noise-free system
for k = 2:T
    % Deterministic part
    fvec = epileptor(taux(k-1),x(:,k-1),x0,y0,Irest1,Irest2,tau0,tau2,gamma); 
    taux(k) = taux(k-1) + dt;
    
    % First subsystem (x1,y1)
    x(1,k) = x(1,k-1) + fvec(1) * dt;
    x(2,k) = x(2,k-1) + fvec(2) * dt;
    
    % Slow variable z
    x(3,k) = x(3,k-1) + fvec(3) * dt;
    
    % Second subsystem (x2,y2)
    x(4,k) = x(4,k-1) + fvec(4) * dt;
    x(5,k) = x(5,k-1) + fvec(5) * dt;
    
    % Dummy variable g
    x(6,k) = x(6,k-1) + fvec(6) * dt;
end
x = x';

% Down-sampling
t = taux(1:round(st/dt):end);
s = x(1:round(st/dt):end,:);

% Measured time-series data
y = s(:,1) + s(:,4);    % x_1 + x_4
y_1 = s(:,1);           % x_1
y_2 = s(:,4);           % x_4
y_3 = s(:,3);           % x_3

% Noise-free time-series data
z = s(:,3);

% Noise
switch noise_type
    case 1
        eta = cos(10000*t) - 0.5*sin(20000*t) + 2*cos(70000*t);
    case 2
        eta_1 = cos(10000*t + 0.7791);
        eta_2 = zeros(length(t),1);
        for i=1:length(t)
            eta_2(i) = max(-100,min(100,0.0375*(sin(100*t(i)))^2*abs(cos(100*t(i))^(-1/2))*sign(cos(100*t(i)) -0.075*(cos(100*t(i)))^2*abs(cos(100*t(i))^(1/2))*sign(cos(100*t(i))))));
        end
        eta = eta_1 + eta_2';
    case 3
        noise_var = 0.01;
        rng(3,'twister')   % for repeatability
        eta = wgn(1,length(t),noise_var,'linear');
    case 4
        noise_var = 0.1;
        rng(4,'twister')   % for repeatability
        eta = wgn(1,length(t),noise_var,'linear');
    case 5
        noise_var = 1;
        rng(5,'twister')   % for repeatability
        eta = wgn(1,length(t),noise_var,'linear');
    otherwise
        disp('Noise type not valid.');
        return;
end

% Measured time-series data
switch additive_noise
    case 0
        y_noisy = y .*(eta'+ones(size(y,1),size(y,2)));
        y_1_noisy = y_1 .*(eta'+ones(size(y,1),size(y,2)));
        y_2_noisy = y_2 .*(eta'+ones(size(y,1),size(y,2)));
        y_3_noisy = y_3 .*(eta'+ones(size(y,1),size(y,2)));
    case 1
        y_noisy = y + eta';
        y_1_noisy = y_1 + eta';
        y_2_noisy = y_2 + eta';
        y_3_noisy = y_3 + eta';
    otherwise
        disp('Data corruption mode not valid.');
        return;
end

signal = y_noisy;

% Savitzky-Golay filtering of Differentiator output
m = 2;
fl = 4551;
fl2 = 2751;

% Estimation of base signal and successive derivatives
t_start_D = tic;
y_D = differentiator(y_noisy,t, n_d, n_f, L_0, q,0);
y_1_D = differentiator(y_1_noisy,t, n_d, n_f, L_0, q,0);
y_2_D = differentiator(y_2_noisy,t, n_d, n_f, L_0, q,0);
y_3_D = differentiator(y_3_noisy,t, n_d, n_f, L_0, q,0);
t_end_D = toc(t_start_D);

y_D_f = sgolayfilt(y_D,m,fl2);
y_1_D_f = sgolayfilt(y_1_D(:,1),m,fl);
y_2_D_f = sgolayfilt(y_2_D(:,1),m,fl);
y_3_D_f = sgolayfilt(y_3_D(:,1),m,fl);
t_end_D_f = toc(t_start_D);

disp(['Computation time Differentiator: ',num2str(round(t_end_D,2,'significant')),' [s]']);
disp(['Computation time Differentiator plus Savitzky-Golay filter: ',num2str(round(t_end_D_f,2,'significant')),' [s]']);

% Elimination Differentiator transitory
if trans_diff > 0
    t = t(round(trans_diff/dt):end);
    s = s(round(trans_diff/dt):end,:);
    y = y(round(trans_diff/dt):end);
    z = z(round(trans_diff/dt):end);
    y_D = y_D(round(trans_diff/dt):end,:);
    y_D_f = y_D_f(round(trans_diff/dt):end,:);
    y_1_D = y_1_D(round(trans_diff/dt):end,:);
    y_2_D = y_2_D(round(trans_diff/dt):end,:);
    y_3_D = y_3_D(round(trans_diff/dt):end,:);
    y_1_D_f = y_1_D_f(round(trans_diff/dt):end,:);
    y_2_D_f = y_2_D_f(round(trans_diff/dt):end,:);
    y_3_D_f = y_3_D_f(round(trans_diff/dt):end,:);
    y_1 = y_1(round(trans_diff/dt):end);
    y_2 = y_2(round(trans_diff/dt):end);
    y_3 = y_3(round(trans_diff/dt):end);
    y_noisy = y_noisy(round(trans_diff/dt):end);
    y_1_noisy = y_1_noisy(round(trans_diff/dt):end);
    y_2_noisy = y_2_noisy(round(trans_diff/dt):end);
    y_3_noisy = y_3_noisy(round(trans_diff/dt):end);
end

% True derivative of y = x_1 + x_4
y_dot = zeros(size(y,1),1);
for i=1:length(y_dot)
    if s(i,1) < 0
        f_1 = s(i,1)^3 - 3 * s(i,1)^2;
    else
        f_1 = (s(i,4) - 0.6 * (s(i,3) - 4)^2) * s(i,1);
    end
    y_dot(i) = s(i,2)-f_1-s(i,3)+Irest1-s(i,5)+s(i,4)-s(i,4)^3+Irest2 +2*s(i,6)-0.3*(s(i,3)-3.5);
end 

%% Plots 
switch noise_type
    case 1
        color_vec = [0.4660 0.6740 0.1880];
    case 2
        color_vec = [0.6350 0.0780 0.1840];
    case 3
        color_vec = [0 0.4470 0.7410];
    case 4
        color_vec = [0.8500 0.3250 0.0980];
    case 5
        color_vec = [0.4940 0.1840 0.5560];
end

% Measured and noise-free signal time-series of x_1
figID = figID + 1;
figure(figID)
plot(t,y_1_noisy,'LineWidth',1,'Color',[color_vec(1),color_vec(2),color_vec(3),0.7]);
hold on;
plot(t,y_1,'k','LineWidth',1.5);
xlim([t(1),t(end)]);
ylabel('$x_1$','interpreter','latex')
legend('Measured signal','Noise-free signal');
ax = gca;
ax.FontSize = 35;
pbaspect([1.2,1,1])

% Measured and noise-free signal time-series of x_4
figID = figID + 1;
figure(figID)
plot(t,y_2_noisy,'LineWidth',1,'Color',[color_vec(1),color_vec(2),color_vec(3),0.7]);
hold on;
plot(t,y_2,'k','LineWidth',1.5);
xlim([t(1),t(end)]);
ylabel('$x_4$','interpreter','latex')
legend('Measured signal','Noise-free signal');
ax = gca;
ax.FontSize = 35;
pbaspect([1.2,1,1])

% Measured and noise-free signal time-series of x_3
figID = figID + 1;
figure(figID)
plot(t,y_3_noisy,'LineWidth',1,'Color',[color_vec(1),color_vec(2),color_vec(3),0.7]);
hold on;
plot(t,y_3,'k','LineWidth',1.5);
xlim([t(1),t(end)]);
ylabel('$x_3$','interpreter','latex')
legend('Measured signal','Noise-free signal');
ax = gca;
ax.FontSize = 35;
pbaspect([1.2,1,1])

% Noise-free and Differentiator output for x_1
figID = figID + 1;
figure(figID)
plt = plot(t,y_1,'LineWidth',1.5);
hold on;
scatter(t,y_1_D(:,1),0.5);
lgd = legend('$x_1$','$x_{1,D}$');
set(lgd, 'interpreter', 'latex','LineWidth',1.5)
uistack(plt,'top');
ax = gca;
ax.FontSize = 35;
xlim([t(1),t(end)])

% Noise-free and Differentiator output for x_4
figID = figID + 1;
figure(figID)
plt = plot(t,y_2,'LineWidth',1.5);
hold on;
scatter(t,y_2_D(:,1),0.5);
xlim([t(1),t(end)])
lgd = legend('$x_4$','$x_{4,D}$');
set(lgd, 'interpreter', 'latex','LineWidth',1.5)
uistack(plt,'top');
ax = gca;
ax.FontSize = 35;

% Noise-free and Differentiator output for x_3
figID = figID + 1;
figure(figID)
plt = plot(t,y_3,'LineWidth',1.5);
hold on;
scatter(t,y_3_D(:,1),0.5);
xlim([t(1),t(end)])
lgd = legend('$x_3$','$x_{3,D}$');
set(lgd, 'interpreter', 'latex','LineWidth',1.5)
uistack(plt,'top');
ax = gca;
ax.FontSize = 35;

% Original state space and Differentiator reconstruction
figID = figID + 1;
figure(figID)
subplot(1,3,1)
scatter3(y_1,y_2,y_3,0.7,[0,0,0]);
xlabel('$x_1$','interpreter','latex'); 
ylabel('$x_4$','interpreter','latex'); 
zlabel('$x_3$','interpreter','latex');
title(['Noise-free',newline,'dynamics']);
view([-33.9,11.4])
ax = gca;
ax.FontSize = 25;
subplot(1,3,2)
scatter3(y_1_noisy,y_2_noisy,y_3_noisy,0.5,color_vec);
xlabel('$y_1$','interpreter','latex'); 
ylabel('$y_2$','interpreter','latex'); 
zlabel('$y_3$','interpreter','latex');
title(['Measured',newline,'dynamics']);
view([-33.9,11.4])
ax = gca;
ax.FontSize = 25;
subplot(1,3,3)
scatter3(y_1_D(:,1),y_2_D(:,1),y_3_D(:,1),0.5,color_vec);
xlabel('$x_{1,Df}$','interpreter','latex'); 
ylabel('$x_{1,Df}$','interpreter','latex'); 
zlabel('$x_{1,Df}$','interpreter','latex');
title(['Differentiator plus SG filter',newline,'reconstruction']);
view([-33.9,11.4])
ax = gca;
ax.FontSize = 25;

% Relative errors
figID = error_comparison_HR_and_Epileptor([y_1,y_2,y_3],[y_1_D(:,1),y_2_D(:,1),y_3_D(:,1)],[y_1_D_f,y_2_D_f,y_3_D_f],40,figID);

% 2D differential embedding
figID = figID + 1;
figure(figID);
subplot(1,3,1); scatter(y,y_dot,0.7,[0,0,0]);
xlabel('$y$','interpreter','latex'); 
ylabel('$\dot{y}$','interpreter','latex');
title(['Noise-free',newline,'dynamics']);
ax = gca;
ax.FontSize = 25;
subplot(1,3,2); scatter(y_D(:,1),y_D(:,2),0.7,color_vec);
xlabel('$y_D$','interpreter','latex'); 
ylabel('$\dot{y}_D$','interpreter','latex');
title(['Differentiator',newline,'reconstruction']);
ax = gca;
ax.FontSize = 25;
subplot(1,3,3); scatter(y_D_f(:,1),y_D_f(:,2),0.7,color_vec);
xlabel('$y_{D,f}$','interpreter','latex'); 
ylabel('$\dot{y}_{D,f}$','interpreter','latex');
title(['Differentiator plus SG filter',newline,'reconstruction']);
ax = gca;
ax.FontSize = 25;

