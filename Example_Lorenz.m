% Uros Sutulovic, 03/2025

clear; close all; clc;
addpath([pwd,'/Systems/']) 
addpath([pwd,'/Util/'])
figID = 0;

noise_type = 2;             % 1 = harmonic
                            % 2 = unbounded
                            % 3 = white Gaussian variance 0.01
                            % 4 = white Gaussian variance 0.1
                            % 5 = white Gaussian variance 1

additive_noise = 1;         % 0 = multiplicative noise
                            % 1 = additive noise

include_grassberger = 0;    % 0 = Schreiber-Grassberger method is skipped
                            % 1 = Schreiber-Grassberger method is included

% Differentiator parameters
n_d = 2;
n_f = 3;            
L_0 = 3.75e4;       
q = 0;              % NaN to disable discrete L-adaptation
trans_diff = 1;     % differentiator transient to discard 

% Schreiber-Grassberger de-noising parameters
iterations = 1;

% Simulation time
dt = 1e-4;
t0 = 0; 
tf = 1020;
tspan = t0:dt:tf;
trans = 1000;                                % transiente time
st = 1e-4;                                   % sampling time
N = length(trans/dt:round(st/dt):tf/dt);     % data points

% Delay embedding parameters
dE = 6;
tau = 1000;

%% Lorenz system simulation
% Parameters
sigma = 10; rho = 28; beta = 8/3;

% Numerical simulation
s0 = [1;1;1];
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
[taux,saux] = ode45(@(t,s)lorenz(t,s,sigma,rho,beta),tspan,s0,options);

% Remove transient, samples data
t = taux(round(trans/dt):round(st/dt):round(tf/dt));
s = saux(round(trans/dt):round(st/dt):round(tf/dt),:);

%% Differentiator-based signal estimation
% Noise-free time-series data
y = s(:,1);

% Data corruption
switch noise_type
    case 1
        eta = cos(10000*t) - 0.5*sin(20000*t) + 2*cos(70000*t);
        noise_var = var(eta);
    case 2
        eta_1 = cos(10000*t + 0.7791);
        eta_2 = zeros(length(t),1);
        for i=1:length(t)
            eta_2(i) = max(-100,min(100,0.0375*(sin(100*t(i)))^2*abs(cos(100*t(i))^(-1/2))*sign(cos(100*t(i)) -0.075*(cos(100*t(i)))^2*abs(cos(100*t(i))^(1/2))*sign(cos(100*t(i))))));
        end
        eta = eta_1 + eta_2;
        noise_var = var(eta);
    case 3
        noise_var = 0.01;
        rng(3,'twister')   % for repeatability
        eta = wgn(length(y),1,noise_var,'linear');
    case 4
        noise_var = 0.1;
        rng(4,'twister')   % for repeatability
        eta = wgn(length(y),1,noise_var,'linear');
    case 5
        noise_var = 1;
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
switch noise_type
    case 1
        fl = 551;
    case 2
        fl = 1999;
    case 3
        fl = 1999;
    case 4
        fl = 2501;
    case 5
        fl = 3001;
end

% Estimation of base signal and successive derivatives
t_start_D = tic;
z = differentiator(signal,t,n_d,n_f,L_0,q,0);
t_end_D = toc(t_start_D);
z_f = sgolayfilt(z(:,1:2),m,fl);
t_end_D_f = toc(t_start_D);

disp(['Computation time Differentiator: ',num2str(round(t_end_D,2,'significant')),' [s]']);
disp(['Computation time Differentiator plus Savitzky-Golay filter: ',num2str(round(t_end_D_f,2,'significant')),' [s]']);

% Elimination differentiator transitory
if trans_diff > 0
    z = z(round(trans_diff/dt):end,:);
    z_f = z_f(round(trans_diff/dt):end,:);
    s = s(round(trans_diff/dt):end,:);
    t = t(round(trans_diff/dt):end);
    signal = signal(round(trans_diff/dt):end);
    y = y(round(trans_diff/dt):end);
end

% Grassberger denoising
switch include_grassberger
    case 0
    case 1
    t_start_G = tic;
    y_denoised = Schreiber_Grassberger(signal,iterations,dE,tau,noise_var);
    t_end_G = toc(t_start_G);
    disp(newline);
    disp(['Computation time Schreiber-Grassberger de-noising: ',num2str(round(t_end_G/3600,2,'significant')),' [h]']);
    otherwise
        disp('Option to include Schreiber-Grassberger method not valid.');
        return;
end

%% Time-delay embeddings
% Noise-free signal
emb{1} = y((dE-1)*tau+1:end,1);
for i = 2:dE
    emb{i} = y((dE-i)*tau+1:end-(i-1)*tau,1);
end
y_delayed = [];
for i = 1:dE
    y_delayed(:,i) = emb{i};
end

clear 'emb';

% Measured signal
emb{1} = signal((dE-1)*tau+1:end,1);
for i = 2:dE
    emb{i} = signal((dE-i)*tau+1:end-(i-1)*tau,1);
end
signal_delayed = [];
for i = 1:dE
    signal_delayed(:,i) = emb{i};
end

clear 'emb';

% Differentiator reconstruction
y_D = z(:,1);
emb{1} = y_D((dE-1)*tau+1:end,1);
for i = 2:dE
    emb{i} = y_D((dE-i)*tau+1:end-(i-1)*tau,1);
end
y_D_delayed = [];
for i = 1:dE
    y_D_delayed(:,i) = emb{i};
end

% Differentiator with SG filter reconstruction
y_D = z_f(:,1);
emb{1} = y_D((dE-1)*tau+1:end,1);
for i = 2:dE
    emb{i} = y_D((dE-i)*tau+1:end-(i-1)*tau,1);
end
y_D_f_delayed = [];
for i = 1:dE
    y_D_f_delayed(:,i) = emb{i};
end

%% Plots
% Measured and noise-free signal time-series
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

% Noise-free signal and reconstructed signal with Differentiator
figID = figID + 1;
figure(figID);
plt = plot(t,s(:,1),'LineWidth',2);
hold on;
plot(t,z(:,1),'LineWidth',2);
lgd = legend('$y$','$y_{D}$');
set(lgd, 'interpreter', 'latex','LineWidth',2)
uistack(plt,'top');
xlim([t(1),t(end)]);
ax = gca;
ax.FontSize = 35;

% Noise-free derivative of the signal and reconstructed derivative of the
% signal with Differentiator
figID = figID + 1;
figure(figID);
plt = plot(t,sigma*(s(:,2)-s(:,1)),'LineWidth',2);
hold on;
plot(t,z(:,2),'LineWidth',2);
lgd = legend('$\dot{y}$','$\dot{y}_{D}$');
set(lgd, 'interpreter', 'latex','LineWidth',2)
uistack(plt,'top');
xlim([t(1),t(end)]);
ax = gca;
ax.FontSize = 35;

% Original coordinates noise-free attractor and Differentiator reconstruction
figID = figID + 1;
figure(figID);
subplot(1,3,1); scatter(s(:,1),s(:,2),0.5,[0,0,0]);
xlabel('$x_1$','interpreter','latex'); ylabel('$x_2$','interpreter','latex');
title(['Noise-free',newline,'dynamics']);
ax = gca;
ax.FontSize = 25;
subplot(1,3,2); scatter(z(:,1),z(:,1) + 1/sigma*z(:,2),0.7,color_vec);
xlabel('$x_{1,D}$','interpreter','latex'); ylabel('$x_{2,D}$','interpreter','latex');
title(['Differentiator',newline,'reconstruction']);
ax = gca;
ax.FontSize = 25;
subplot(1,3,3); scatter(z_f(:,1),z_f(:,1) + 1/sigma*z_f(:,2),0.7,color_vec);
xlabel('$x_{1,Df}$','interpreter','latex'); ylabel('$x_{2,Df}$','interpreter','latex');
title(['Differentiator plus SG filter',newline,'reconstruction']);
ax = gca;
ax.FontSize = 25;

% Time-delay embedding noise-free dynamics, measured dynamics and de-noised dynamics
figID = figID + 1;
figure(figID);
switch include_grassberger
    case 0
        subplot(1,3,1); 
    case 1
        subplot(1,4,1); 
end
scatter(y_delayed(:,1),y_delayed(:,2),0.4,[0,0,0]);
xlabel('$x_1(t)$','interpreter','latex'); ylabel('$x_1(t+\tau)$','interpreter','latex');
title(['Takens embedding',newline,'noise-free signal']);
ax = gca;
ax.FontSize = 25;
switch include_grassberger
    case 0
        subplot(1,3,2); 
    case 1
        subplot(1,4,2); 
end 
scatter(signal_delayed(:,1),signal_delayed(:,2),0.4,color_vec);
xlabel('$y(t)$','interpreter','latex'); ylabel('$y(t+\tau)$','interpreter','latex');
title(['Takens embedding',newline,'measured signal']);
ax = gca;
ax.FontSize = 25;
switch include_grassberger
    case 0
        subplot(1,3,3); 
    case 1
        subplot(1,4,3); 
end 
scatter(y_D_delayed(:,1),y_D_delayed(:,2),0.4,color_vec);
xlabel('$y_D(t)$','interpreter','latex'); ylabel('$y_D(t+\tau)$','interpreter','latex');
title(['Takens embedding',newline,'Differentiator']);
ax = gca;
ax.FontSize = 25;
switch include_grassberger
    case 1
        subplot(1,4,4); scatter(y_denoised(:,1),y_denoised(:,2),0.4,color_vec);
        xlabel('$y_G(t)$','interpreter','latex'); ylabel('$y_G(t+\tau_d)$','interpreter','latex');
        title(['Takens embedding',newline,'Schreiber-Grassberger']);
        ax = gca;
        ax.FontSize = 25; 
end

switch include_grassberger
    case 0
        % Relative errors
        figID = error_comparison_Lorenz(y_delayed,y_D_delayed,y_D_f_delayed,NaN,40,figID);
        
        % Space-filling patterns
        figID = space_filling_Lorenz(y_delayed, signal_delayed, y_D_delayed, NaN, [4.8,11.75;4.8,11.75], [-11.46,-5.15;-11.46,-5.15],figID);
    case 1
        % Relative errors
        figID = error_comparison_Lorenz(y_delayed,y_D_delayed,y_D_f_delayed,y_denoised,40,figID);
        
        % Space-filling patterns
        figID = space_filling_Lorenz(y_delayed, signal_delayed, y_D_delayed, y_denoised, [4.8,11.75;4.8,11.75], [-11.46,-5.15;-11.46,-5.15],figID);
end

