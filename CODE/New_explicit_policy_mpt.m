clear
clc

% parameters
g = 9.8; % gravitational coefficient
m = 0.15; % mass
l = 0.5; % length
mu = 0.05; % frictional coefficient
dt = 0.02; % sampling period

AG = [1,      dt;...
    g/l*dt, 1-mu/(m*l^2)*dt];
BG = [0; dt/(m*l^2)];
bar_varphi=0.73;
nG = size(AG, 1);
nu = size(BG, 2);
ls = (bar_varphi-sin(bar_varphi))/bar_varphi;
ms = 0; 
[A,B,C,D,statenum,colnum_A_Phi,omeganum,M_Theta] = Pendulum(dt,m,l,mu,g,ls,ms); 
% Linear discrete-time prediction model
model=LTISystem('A', AG, 'B', BG);

% Input constraints
model.u.min =-2;
model.u.max = 2;

% State constraints
model.x.min = [-2.5; -6];
model.x.max = [2.5; 6];

% constraint sets represented as polyhedra
X = Polyhedron('lb',model.x.min,'ub',model.x.max);
U = Polyhedron('lb',model.u.min,'ub',model.u.max);

% Penalties in the cost function
Q = 10e8*eye(nG);
R = 10e-8;
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);

% Maximal Invariant Set Computation
[Pinf,Kinf,L] = idare(AG,BG,Q,R);% closed loop system
Acl=AG-BG*Kinf;
S=X.intersect(Polyhedron('H',[-U.H(:,1:nu)*Kinf U.H(:,nu+1)]));
Oinf=max_pos_inv(Acl,S);

model.x.with('terminalSet');
model.x.terminalSet = Oinf;
model.x.with('terminalPenalty');
model.x.terminalPenalty = QuadFunction(Pinf);
umax=0.67;
% Online MPC object
online_ctrl = MPCController( model, 3 );


explicit_ctrl = online_ctrl.toExplicit();


if false
    hold on
    
    num_pt = 100;
    xU = model.x.max;
    xL = model.x.min;
    x1 = linspace(xL(1), xU(1),num_pt);
    x2 = linspace(xL(2),xU(2),num_pt);
    ug = zeros(num_pt,num_pt);
    data = [];
    for i=1:num_pt
        for j = 1:num_pt
            ug(i, j) = online_ctrl.evaluate([x1(i); x2(j)]);
            if ~isnan(ug(i, j))
                data = [data; x1(i), x2(j), ug(i, j)];
            end
        end
    end
    mesh(x1,x2,ug')
    
    % save the state/control data pairs
    writematrix(data,'exp_data.csv')
end

%% simulate
% initial state

if true
    x0 =  [0.23,3.5];
    %[0.23,3.5]; [0.43,3.0];  [-0.33, -3.3];
    eta0 = [x0.';0];
    etaa{1} = eta0;
    T = 16; %second
    x = x0.'; 
    time = 0:dt:T;
    X_traj = zeros(2, length(time));
    U_traj = zeros(1, length(time));
    xi_traj = zeros(1, length(time));
    process_noise_std = 2*10e-4;
    for t = 1:length(time)
        X_traj(:, t) = x; 
        u = explicit_ctrl.evaluate(x); 
        u = max(min(u, umax), -umax);
        U_traj(t) = u; 
        etaa{t}(1)=x(1);
        etaa{t}(2)=x(2);
        etaa{t+1} = A*etaa{t}+B*[u; x(1)-sin(x(1))];
        x = AG * x + BG * u;  % 加入噪声项
        xi_traj(t) = etaa{t}(3);
    end
    folderName = sprintf('%.2f_%.2f', x0(1), x0(2));
    folderName = strrep(folderName, '-', 'm');
    
    if ~isfolder(folderName)
        mkdir(folderName); % 如果不存在，则创建
    end

    writematrix(X_traj,fullfile(folderName, 'expert_state.csv'))
    writematrix(xi_traj,fullfile(folderName, 'expert_xi.csv'))
    
    figure;
    subplot(2, 1, 1);
    plot(time, X_traj(1, :), 'LineWidth', 2);
    hold on;
    title('state_traj');
    xlabel('Sec.');
    ylabel('state');
    legend(' x1');
    
    subplot(2, 1, 2);
    plot(time, U_traj, 'LineWidth', 2);
    title('U_traj');
    xlabel('Sec.');
    ylabel('input');
end