clear
clc

% parameters
g = 9.8; % gravitational coefficient
m = 0.15; % mass
l = 0.5; % length
mu = 0.05; % frictional coefficient
dt = 0.05; % sampling period (adjusted to make convergence slower)

AG = [1,      dt;...
    g/l*dt, 1-mu/(m*l^2)*dt];
BG = [0; dt/(m*l^2)];
bar_varphi = 0.73;
nG = size(AG, 1);
nu = size(BG, 2);
ls = (bar_varphi - sin(bar_varphi)) / bar_varphi;
ms = 0;
[A,B,C,D,statenum,colnum_A_Phi,omeganum,M_Theta] = Pendulum(dt, m, l, mu, g, ls, ms); 
% Linear discrete-time prediction model
model = LTISystem('A', AG, 'B', BG);

% Input constraints
model.u.min = -0.7; % 控制输入下限
model.u.max = 0.7;  % 控制输入上限

% State constraints
model.x.min = [-2.5; -6];
model.x.max = [2.5; 6];

% constraint sets represented as polyhedra
X = Polyhedron('lb',model.x.min,'ub',model.x.max);
U = Polyhedron('lb',model.u.min,'ub',model.u.max);

% Penalties in the cost function (adjusted for longer convergence)
Q = 2 * eye(nG); % Increased state penalty
R = 1; % Decreased control input penalty
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);

% Maximal Invariant Set Computation
[Pinf,Kinf,L] = idare(AG,BG,Q,R); % closed loop system
Acl = AG - BG * Kinf;
S = X.intersect(Polyhedron('H',[-U.H(:,1:nu) * Kinf U.H(:,nu+1)]));
Oinf = max_pos_inv(Acl, S);

model.x.with('terminalSet');
model.x.terminalSet = Oinf;
model.x.with('terminalPenalty');
model.x.terminalPenalty = QuadFunction(Pinf);

% Online MPC object (with increased prediction horizon)
online_ctrl = MPCController(model, 10); % Increased prediction horizon

explicit_ctrl = online_ctrl.toExplicit();

% simulate
% initial state
x0 = [0.23, 3.5];
eta0 = [x0.';0];
etaa{1} = eta0;
T = 16; % second
x = x0.'; 
time = 0:dt:T;

X_traj = zeros(2, length(time));
U_traj = zeros(1, length(time));
xi_traj = zeros(1, length(time));
for t = 1:length(time)
    % 加入测量噪声
    noisy_x = x + 0.05 * randn(size(x));  % 示例：加入测量噪声
    
    X_traj(:, t) = noisy_x;  % 使用带噪声的观测值
    u = explicit_ctrl.evaluate(noisy_x);  % 根据带噪声的观测值计算控制输入
    U_traj(t) = u; 
    
    % 保证 u 在[-0.7, 0.7]范围内
    u = max(min(u, 0.7), -0.7);  % 将控制输入限制在 [-0.7, 0.7] 范围内

    etaa{t}(1) = x(1);
    etaa{t}(2) = x(2);
    etaa{t+1} = A * etaa{t} + B * [u; x(1) - sin(x(1))];

    % 系统状态更新
    x = AG * x + BG * u; 
    
    xi_traj(t) = etaa{t}(3);
end

folderName = sprintf('%.2f_%.2f', x0(1), x0(2));
folderName = strrep(folderName, '-', 'm');

if ~isfolder(folderName)
    mkdir(folderName);
end

writematrix(X_traj, fullfile(folderName, 'expert_state.csv'));
writematrix(xi_traj, fullfile(folderName, 'expert_xi.csv'));

figure;
subplot(2, 1, 1);
plot(time, X_traj(1, :), 'LineWidth', 2);
hold on;
plot(time, X_traj(2, :), 'LineWidth', 2);
title('state_traj');
xlabel('Sec.');
ylabel('state');
legend('x1', 'x2');

subplot(2, 1, 2);
plot(time, xi_traj, 'LineWidth', 2);
title('xi_traj');
xlabel('Sec.');
ylabel('input');
