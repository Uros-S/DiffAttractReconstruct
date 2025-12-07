% Uros Sutulovic, 03/2025
% If you use this code, please reference:
% U Sutulovic, D Proverbio, R Katz, G Giordano. 
% Efficient and faithful reconstruction of dynamical attractors using homogeneous differentiators. 
% Chaos, Solitons & Fractals 199 (3), 116798 (2025)

clear; close all; clc;
addpath([pwd,'/Systems/'])
addpath([pwd,'/Util/'])
figID = 0;

regime = 4;                         % 1 = quiescence 
                                    % 2 = spiking
                                    % 3 = square-wave
                                    % 4 = plateau
                                    % 5 = chaotic

noise_type = 3;                     % 1 = harmonic
                                    % 2 = unbounded
                                    % 3 = white Gaussian variance 0.001
                                    % 4 = white Gaussian variance 0.01
                                    % 5 = white Gaussian variance 0.1

additive_noise = 0;                 % 0 = multiplicative noise
                                    % 1 = additive noise

% Differentiator parameters
n_d = 3;
n_f = 9;
L_0 = 10;
q = 0;              % NaN to disable discrete L-adaptation
trans_diff = 5;     % differentiator transient to discard

% Simulation time
dt = 1e-4; 
t0 = 0; 
tf = 2500;  
tspan = t0:dt:tf;
trans = 1500;                                 % transiente time
st = 1e-4;                                    % sampling time
N = length(trans/dt:round(st/dt):tf/dt);      % data points

 %% Hindmarsh-Rose system simulation  
% Parameters
a = 1; c = 1; d = 5; r = 0.01; sdyn = 4; xR = -(1+sqrt(5))/2;
switch regime
    case 1 
        I = 2.2; b = 3.2; 
    case 2
        I = 2.5; b = 3;
    case 3
        I = 3; b = 2.7;
    case 4
        I = 4; b = 2.5;
    case 5
        I = 2.70; b = 3.03;
    otherwise
        disp('Dynamical regime selected not valid.');
        return;
end

% Numerical simulation
s0 = [0.65; 0.55; 0.45];
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
[taux,saux] = ode45(@(t,s)neuronHR(t,s,a,b,c,d,I,r,sdyn,xR),tspan,s0,options);

% Remove transient, samples data
t = taux(round(trans/dt):round(st/dt):round(tf/dt));
s = saux(round(trans/dt):round(st/dt):round(tf/dt),:);

%% Differentiator-based signal estimation
% Noise-free time-series data
y = s(:,3);

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
        eta = eta_1 + eta_2;
    case 3
        noise_var = 0.001;
        rng(3,'twister')   % for repeatability
        eta = wgn(length(y),1,noise_var,'linear');
    case 4
        noise_var = 0.01;
        rng(4,'twister')   % for repeatability
        eta = wgn(length(y),1,noise_var,'linear');
    case 5
        noise_var = 0.1;
        rng(5,'twister')   % for repeatability
        eta = wgn(length(y),1,noise_var,'linear');
    otherwise
        disp('Noise type not valid.');
        return;
end

% Measured time-series data
switch additive_noise
    case 0
        signal = y.*(eta+ones(size(eta,1),size(eta,2)));
    case 1
        signal = y + eta;
    otherwise
        disp('Data corruption mode not valid.');
        return;
end

% Savitzky-Golay filtering of Differentiator output
m = 2;
fl = 8251;

% Estimation of base signal and successive derivatives
t_start_D = tic;
z = differentiator(signal,t,n_d,n_f,L_0,q,0);
t_end_D = toc(t_start_D);
z_f = sgolayfilt(z(:,1:3),m,fl);
t_end_D_f = toc(t_start_D);

disp(['Computation time Differentiator: ',num2str(round(t_end_D,2,'significant')),' [s]']);
disp(['Computation time Differentiator plus Savitzky-Golay filter: ',num2str(round(t_end_D_f,2,'significant')),' [s]']);

% Elimination Differentiator transitory
if trans_diff > 0
    z = z(round(trans_diff/dt):end,:);
    z_f = z_f(round(trans_diff/dt):end,:);
    s = s(round(trans_diff/dt):end,:);
    t = t(round(trans_diff/dt):end);
    signal = signal(round(trans_diff/dt):end);
    y = y(round(trans_diff/dt):end);
end

% State reconstruction
x_1_D = xR*ones(length(z(:,1)),1) + 1/sdyn*(z(:,1) + 1/r*(z(:,2)));
x_2_D = 1/sdyn*(z(:,2) + 1/r*z(:,3)) + a*x_1_D.^3 - b*x_1_D.^2 + z(:,1) - I*ones(length(z(:,1)),1);
x_3_D = z(:,1);

x_1_D_f = xR*ones(length(z_f(:,1)),1) + 1/sdyn*(z_f(:,1) + 1/r*(z_f(:,2)));
x_2_D_f = 1/sdyn*(z_f(:,2) + 1/r*z_f(:,3)) + a*x_1_D_f.^3 - b*x_1_D_f.^2 + z_f(:,1) - I*ones(length(z_f(:,1)),1);
x_3_D_f = z_f(:,1);

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

% Measured and noise-free signal time-series
figID = figID + 1;
figure(figID);
plot(t,signal,'LineWidth',1,'Color',[color_vec(1),color_vec(2),color_vec(3),0.7]);
hold on;
plot(t,y,'k','LineWidth',1.5);
legend('Measured signal','Noise-free signal');
xlim([t(1),t(end)]);
pbaspect([2,1,1])
ax = gca;
ax.FontSize = 35;

% Noise-free signal and reconstructed signal
figID = figID + 1;
figure(figID);
plt = plot(t,y,'LineWidth',1.5);
hold on;
plot(t,z(:,1),'LineWidth',1.5);
lgd = legend('$y$','$y_{D}$');
set(lgd, 'interpreter', 'latex','LineWidth',1.5)
uistack(plt,'top');
ax = gca;
ax.FontSize = 35;

% Noise-free first derivative of the signal and reconstructed first 
% derivative of the signal with Differentiator
figID = figID + 1;
figure(figID);
plt = plot(t,r*(sdyn*(s(:,1)-xR*ones(length(s(:,1)),1))-s(:,3)),'LineWidth',1.5);
hold on;
plot(t,z(:,2),'LineWidth',1.5);
lgd = legend('$\dot{y}$','$\dot{y}_D$');
set(lgd, 'interpreter', 'latex','LineWidth',1.5)
uistack(plt,'top');
ax = gca;
ax.FontSize = 35;

% Noise-free second derivative of the signal and reconstructed second 
% derivative of the signal with Differentiator
figID = figID + 1;
figure(figID);
plt = plot(t,r*sdyn*(s(:,2) - a*s(:,1).^3 + b*s(:,1).^2 - s(:,3) + I*ones(length(s(:,1)),1)) - ...
       r^2*(sdyn*(s(:,1) - xR*ones(length(s(:,1)),1)) - s(:,3)),'LineWidth',1.5);
hold on;
plot(t,z(:,3),'LineWidth',1.5);
lgd = legend('$\ddot{y}$','$\ddot{y}_D$');
set(lgd, 'interpreter', 'latex','LineWidth',1.5)
uistack(plt,'top');
ax = gca;
ax.FontSize = 35;

% Original state space and Differentiator reconstruction
figID = figID + 1;
figure(figID);
subplot(1,3,1); scatter3(s(:,1),s(:,2),s(:,3),0.5,[0,0,0]);
xlabel('$x_1$','interpreter','latex'); 
ylabel('$x_2$','interpreter','latex'); 
zlabel('$x_3$','interpreter','latex');
title(['Noise-free',newline,'dynamics']);
view([-14,17])
ax = gca;
ax.FontSize = 25;
subplot(1,3,2); scatter3(x_1_D,x_2_D,x_3_D,0.7,color_vec);
xlabel('$x_{1,D}$','interpreter','latex'); 
ylabel('$x_{2,D}$','interpreter','latex'); 
zlabel('$x_{3,D}$','interpreter','latex');
title(['Differentiator',newline,'reconstruction']);
view([-14,17])
ax = gca;
ax.FontSize = 25;
subplot(1,3,3); scatter3(x_1_D_f,x_2_D_f,x_3_D_f,0.7,color_vec);
xlabel('$x_{1,Df}$','interpreter','latex'); 
ylabel('$x_{2,Df}$','interpreter','latex'); 
zlabel('$x_{3,Df}$','interpreter','latex');
title(['Differentiator plus SG filter',newline,'reconstruction']);
view([-14,17])
ax = gca;
ax.FontSize = 25;

% Relative errors
figID = error_comparison_HR([s(:,1),s(:,2),s(:,3)],[x_1_D,x_2_D,x_3_D],[x_1_D_f,x_2_D_f,x_3_D_f],40,figID);

