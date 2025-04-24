clear
clc

% Initial Parameters
Nstep = 800;                % Number of time steps
Ts = 0.02;                  % Sampling time
 
% Time vector
time = Ts * (0:Nstep-1); 
%% Initial values case 1
x_0c1 = [0.23, 3.5];        
folder = sprintf('%.2f_%.2f/', x_0c1(1), x_0c1(2)); 
folder = strrep(folder, '-', 'm');  % Replace '-' with 'minus'

% Load data from CSV files
% ET

etaa1_c1 = load([folder 'etaa41.csv']);
etaa2_c1 = load([folder 'etaa42.csv']); 
etaa3_c1 = load([folder 'etaa43.csv']); 

%ST
staa1_c1 = load([folder 'staa1.csv']);
staa2_c1 = load([folder 'staa2.csv']);
staa3_c1 = load([folder 'staa3.csv']);

%Expert
exp1 = load([folder 'expert_state.csv']);
exp_xi1 = load([folder 'expert_xi.csv']);
%% Initial values case 2
x_0c2 = [0.43, 3.0]; 

folder = sprintf('%.2f_%.2f/', x_0c2(1), x_0c2(2)); 
folder = strrep(folder, '-', 'm');  % Replace '-' with 'minus

% Load data from CSV files
% ET

etaa1_c2 = load([folder 'etaa41.csv']);
etaa2_c2 = load([folder 'etaa42.csv']); 
etaa3_c2 = load([folder 'etaa43.csv']); 

%ST
staa1_c2 = load([folder 'staa1.csv']);
staa2_c2 = load([folder 'staa2.csv']);
staa3_c2 = load([folder 'staa3.csv']);

%Expert
exp2 = load([folder 'expert_state.csv']);
exp_xi2 = load([folder 'expert_xi.csv']);


%% Initial values case 3

x_0c3 = [-0.33, -3.3]; 

folder = sprintf('%.2f_%.2f/', x_0c3(1), x_0c3(2)); 
folder = strrep(folder, '-', 'm');  % Replace '-' with 'minus

% Load data from CSV files
% ET

etaa1_c3 = load([folder 'etaa41.csv']);
etaa2_c3 = load([folder 'etaa42.csv']); 
etaa3_c3 = load([folder 'etaa43.csv']); 

%ST
staa1_c3 = load([folder 'staa1.csv']);
staa2_c3 = load([folder 'staa2.csv']);
staa3_c3 = load([folder 'staa3.csv']);

%Expert
exp3 = load([folder 'expert_state.csv']);
exp_xi3 = load([folder 'expert_xi.csv']);

%% plot ET vs ST
figure(111)

subplot(3,1,1)

plot(time, exp1(1, 1:Nstep), '-', 'linewidth', 1.5)
hold on
plot(time, etaa1_c1, '-', 'linewidth', 1.5);
hold on
plot(time, staa1_c1, '-', 'linewidth', 1.5);
hold off

xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('$x_{1}$ (MPC)', '$x_{1}$ (ET,$\vartheta=1$)', '$x_{1}$ (ST)',  'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered ,Self-Triggered, and MPC State Response of Angular Position $x_{1}(\mathrm{rad})$', 'Interpreter','latex')

subplot(3,1,2)
plot(time, exp1(2, 1:Nstep), '-', 'linewidth', 1.5)
hold on
plot(time, etaa2_c1, '-', 'linewidth', 1.5);
hold on
plot(time, staa2_c1, '-', 'linewidth', 1.5);
hold off
xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('$x_{2}$ (MPC)', '$x_{2}$ (ET,$\vartheta=1$)', '$x_{2}$ (ST)', '$x_{2}$ (MPC)', '$x_{2}$(ET) ','$x_{2}$(ST)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered ,Self-Triggered, and MPC State Response of Angular Velocity $x_{2}(\mathrm{rad/s})$', 'Interpreter', 'latex')


subplot(3,1,3)
plot(time,exp_xi1(1, 1:Nstep),'linewidth',1.5);
hold on
plot(time,etaa3_c1,'-','linewidth',1.5);
hold on 
plot(time,staa3_c1,'-.','linewidth',1.5);
hold off
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('$\xi$(MPC)','$\xi$ (ET,$\vartheta=1$)','$\xi$ (ST)','FontSize',10,'Orientation','horizon');
set(hh, 'Interpreter', 'latex');
t='The Event-Triggered, Self-Triggered and MPC State Response of IQC Virtual-Filter State $\xi(\mathrm{rad})$';
title(t, 'Interpreter', 'latex')



figure(112)

subplot(3,1,1)

plot(time, exp2(1, 1:Nstep), '-', 'linewidth', 1.5)
hold on
plot(time, etaa1_c2, '-', 'linewidth', 1.5);
hold on
plot(time, staa1_c2, '-', 'linewidth', 1.5);
hold off

xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('$x_{1}$ (MPC)', '$x_{1}$ (ET,$\vartheta=1$)', '$x_{1}$ (ST)',  'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered ,Self-Triggered, and MPC State Response of Angular Position $x_{1}(\mathrm{rad})$', 'Interpreter','latex')

subplot(3,1,2)
plot(time, exp2(2, 1:Nstep), '-', 'linewidth', 1.5)
hold on
plot(time, etaa2_c2, '-', 'linewidth', 1.5);
hold on
plot(time, staa2_c2, '-', 'linewidth', 1.5);
hold off
xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('$x_{2}$ (MPC)', '$x_{2}$ (ET,$\vartheta=1$)', '$x_{2}$ (ST)', '$x_{2}$ (MPC)', '$x_{2}$(ET) ','$x_{2}$(ST)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered ,Self-Triggered, and MPC State Response of Angular Velocity $x_{2}(\mathrm{rad/s})$', 'Interpreter', 'latex')


subplot(3,1,3)
plot(time,exp_xi2(1, 1:Nstep),'linewidth',1.5);
hold on
plot(time,etaa3_c2,'-','linewidth',1.5);
hold on 
plot(time,staa3_c2,'-.','linewidth',1.5);
hold off
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('$\xi$(MPC)','$\xi$ (ET,$\vartheta=1$)','$\xi$ (ST)','FontSize',10,'Orientation','horizon');
set(hh, 'Interpreter', 'latex');
t='The Event-Triggered, Self-Triggered and MPC State Response of IQC Virtual-Filter State $\xi(\mathrm{rad})$';
title(t, 'Interpreter', 'latex')

figure(113)

subplot(3,1,1)

plot(time, exp3(1, 1:Nstep), '-', 'linewidth', 1.5)
hold on
plot(time, etaa1_c3, '-', 'linewidth', 1.5);
hold on
plot(time, staa1_c3, '-', 'linewidth', 1.5);
hold off

xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('$x_{1}$ (MPC)', '$x_{1}$ (ET,$\vartheta=1$)', '$x_{1}$ (ST)',  'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered ,Self-Triggered, and MPC State Response of Angular Position $x_{1}(\mathrm{rad})$', 'Interpreter','latex')

subplot(3,1,2)
plot(time, exp3(2, 1:Nstep), '-', 'linewidth', 1.5)
hold on
plot(time, etaa2_c3, '-', 'linewidth', 1.5);
hold on
plot(time, staa2_c3, '-', 'linewidth', 1.5);
hold off
xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('$x_{2}$ (MPC)', '$x_{2}$  (ET,$\vartheta=1$)', '$x_{2}$ (ST)', '$x_{2}$ (MPC)', '$x_{2}$(ET) ','$x_{2}$(ST)', 'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered ,Self-Triggered, and MPC State Response of Angular Velocity $x_{2}(\mathrm{rad/s})$', 'Interpreter', 'latex')


subplot(3,1,3)
plot(time,exp_xi3(1, 1:Nstep),'linewidth',1.5);
hold on
plot(time,etaa3_c3,'-','linewidth',1.5);
hold on 
plot(time,staa3_c3,'-.','linewidth',1.5)
hold off
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('$\xi$(MPC)','$\xi$ (ET,$\vartheta=1$)','$\xi$ (ST)','FontSize',10,'Orientation','horizon');
set(hh, 'Interpreter', 'latex');
t='The Event-Triggered, Self-Triggered and MPC State Response of IQC Virtual-Filter State $\xi(\mathrm{rad})$';
title(t, 'Interpreter', 'latex')


%% plot ET in defferent initial states
figure(222)

subplot(3,1,1)

plot(time, etaa1_c1, '-', 'linewidth', 1.5)
hold on
plot(time, etaa1_c2, '-', 'linewidth', 1.5);
hold on
plot(time, staa1_c3, '-', 'linewidth', 1.5);
hold off

xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('Initial State 1', 'Initial State 2', 'Initial State 3',  'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered State Response of Angular Position $x_{1}(\mathrm{rad})$ in Different Initial States($\vartheta=1$)', 'Interpreter','latex')

subplot(3,1,2)
plot(time, etaa2_c1, '-', 'linewidth', 1.5)
hold on
plot(time, etaa2_c2, '-', 'linewidth', 1.5);
hold on
plot(time, staa2_c3, '-', 'linewidth', 1.5);
hold off

xlabel('Time (Sec.)', 'Interpreter', 'latex')
legend('Initial State 1', 'Initial State 2', 'Initial State 3',  'FontSize', 10, 'Orientation', 'horizontal', 'Interpreter', 'latex')
title('The Event-Triggered State Response of Angular Velocity $x_{2}(\mathrm{rad/s})$ in Different Initial States($\vartheta=1$)', 'Interpreter','latex')


subplot(3,1,3)
plot(time,etaa3_c1,'-','linewidth',1.5);
hold on 
plot(time,etaa3_c2,'-.','linewidth',1.5);
hold on 
plot(time,etaa3_c3,'-.','linewidth',1.5);
hold off
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('Initial State 1','Initial State 2','initial State 3','FontSize',10,'Orientation','horizon');
set(hh, 'Interpreter', 'latex');
t='The Event-Triggered State Response of IQC Virtual-Filter State $\xi(\mathrm{rad})$ in Different Initial States($\vartheta=1$)';
title(t, 'Interpreter', 'latex')

%}