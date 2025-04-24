%%%%%%%%% RoA_Comparision
%%%%%%%%% Renjie Ma, Harbin Institute of Technology
%%%%%%%%% Dec 2023

clear;
clc;

%% System Parameters
Ts = 0.01; 
m = 0.15;
l = 0.5;
mu = 0.5;
g = 9.8;
bar_varphi = 0.73; 
ls = (bar_varphi-sin(bar_varphi))/bar_varphi;
ms = 0; 
[A,B,C,D,statenum,colnum_A_Phi,omeganum,M_Theta] = Pendulum(Ts,m,l,mu,g,ls,ms); 

%% Read Pretrained DNN
folder = 'NN_Parameters/'; % Provided by the open-source code [Yin2021Stability]
load([folder 'W1.csv']) % a1 * a0(=statenum)
load([folder 'W2.csv']) % a2 * a1
load([folder 'W3.csv']) % a3(=inputnum) * a2
W ={W1,W2,W3};
a0 = statenum;
a1 = size(W{1},1);
a2 = size(W{2},1);
a3 = size(W{3},1);
aa ={a1,a2}; % m = Col(m_1,m_2)

bias1 = zeros(a1,1);
bias2 = zeros(a2,1);
bias3 = zeros(a3,1);
bias = {bias1,bias2,bias3}; 
umax = 0.7; 
umin = -umax; % saturation

%% Define the initial values

% augmented state eta(0)=eta0
xi0 = zeros(colnum_A_Phi,1);
%x0 = [0.23;-0.88];
x0 = [0.23;3.5]; % 第一次测试：0.23；3.5 ； 第二次测试： 0.43 -3.0 ； 第三次测试： -0.33；-1.8
%x0 = [0.23;3.5];
eta0 = [x0;xi0];
Nstep = 800;

% Augmented state updates of Time-triggered DNN-controlled system 
[eta,u] = T2NNClosedLoop(Nstep,eta0,W,bias,A,B,umax,umin,statenum); % Adjust the saturation via T2NNClosedLoop function
% save the state/control data pairs
writematrix(eta,'t2eta.csv')  
writematrix(u,'t2u.csv')
      
% figure(1); 
% plot(eta(1,1:Nstep),'-. ','linewidth',1.8);
% hold on
% plot(eta(2,1:Nstep),'--','linewidth',1.3);
% hold on
% plot(eta(3,1:Nstep),'-','linewidth',1.3);
% hold off
% xlabel('Time(Sec.)');
% hh = legend('Closed-Loop State $x_{1}(k)$ (Case 1)  ',...
%     'Closed-Loop State: $x_{2}(k)$ (Case 1) ','Virtual Filter State: $\xi(k)$ (Case 1) ','FontSize',11);
% set(hh,'Interpreter','latex')
% 
% figure(2); 
% plot(u(1,1:Nstep),'- ','linewidth',1.8);
% hh = legend('Neural Network Controller $u(k)$ (Case 1)','FontSize',11);
% set(hh,'Interpreter','latex')

n = statenum; % The row number of system state (vector)
psi = size(A,1) - statenum; % The row number of virtual IQC filter state (vector)
varsigma = n+psi; % The row number of augmented state (vector)
ell = 2; % The number of DNN layers (l-layer DNN)
aaa = 0;
for i = 1:ell 
    aaa = aaa+aa{i};
end
a = aaa; % The row number of augmented inputs of activation layers
w = omeganum; % The row number of omega
bar_n = varsigma+a+w+5*n;
m = size(umax,1);

% Define the equilibrium points
x_eq = zeros(n,1);
p1_eq = W{1}*x_eq+bias{1};
m1_eq = tanh(p1_eq);
p2_eq = W{2}*m1_eq+bias{2};
m2_eq = tanh(p2_eq);
p3_eq = W{3}*m2_eq+bias{3};
u_eq = p3_eq;
p_eq = [p1_eq;p2_eq];
m_eq = [m1_eq;m2_eq];

xi_eq = zeros(varsigma-n,1);
eta_eq = [x_eq;xi_eq];
omega_eq = x_eq(1,:)-sin(x_eq(1,:));
nu_eq = x_eq(1,:);

%% IBP & Quadratic Constraints of DNN

Pi_ux = zeros(m,n);
Pi_um = [zeros(m,a-aa{ell}),W{ell+1}];
Pi_u1 = bias{ell+1};
Pi_px = [W{1};zeros(a-aa{1},n)]; 
tilde_Pi_pm = [W{2},zeros(aa{2},aa{2})];
Pi_pm = [zeros(aa{1},a);tilde_Pi_pm];
Pi_p1 = [bias{1};bias{2}];

count=5; 

Delta=[0.25,0.35,0.45,0.25,0.35,0.45,0.25,0.35,0.45;...
    0.97,0.97,0.97,1,1,1,1.03,1.03,1.03];

delta_rho = Delta(1,count); % This value can be adjusted for comparsion 
delta_sigma = Delta(2,count); % This value can be adjusted for comparsion


p1_upp = p1_eq + delta_rho*ones(a1,1); % The upper bound of p1
p1_low = p1_eq - delta_rho*ones(a1,1); % The lower bound of p1

m1_upp = tanh(p1_upp); % The upper bound of m1
m1_low = tanh(p1_low); % The lower bound of m1

p2_upp = [];
p2_low = [];
cc1 = 0.5*(m1_upp+m1_low); % According to [Eq.(19),Yin2021Stability]
rr1 = 0.5*(m1_upp-m1_low);
for i = 1:a2
    p2_upp_i = W{2}(i,:)*cc1+bias{2}(i,:)+abs( W{2}(i,:)*rr1);
    p2_low_i = W{2}(i,:)*cc1+bias{2}(i,:)-abs( W{2}(i,:)*rr1);
    p2_upp = [p2_upp;p2_upp_i];
    p2_low = [p2_low;p2_low_i];
end

rho_1 = min( (tanh(p1_upp) - tanh(p1_eq))./(p1_upp-p1_eq),...
              (tanh(p1_eq) - tanh(p1_low))./(p1_eq-p1_low)  );
sigma= delta_sigma;

rho_2 = min( (tanh(p2_upp) - tanh(p2_eq))./(p2_upp-p2_eq),...
              (tanh(p2_eq) - tanh(p2_low))./(p2_eq-p2_low)  );

M_rho = blkdiag(diag(rho_1),diag(rho_2));
M_sigma = delta_sigma*eye(a);

mathcalM_DNN1 = [M_sigma, -eye(a);-M_rho, eye(a)]; % This aims at formulate the transpose term embedded in Lemma 2

%% Linear Matrix Inequality (c.t.: Event-Triggered Case(Case 2))

L_1 = [eye(varsigma),zeros(varsigma,bar_n-varsigma)]; 
L_2 = [zeros(a,varsigma),eye(a),zeros(a,bar_n-varsigma-a)]; 
L_3 = [zeros(w,varsigma+a),eye(w),zeros(w,5*n)];
L_4 = [zeros(n,varsigma+a+w),eye(n),zeros(n,4*n)];
L_5 = [zeros(n,varsigma+a+w+n),eye(n),zeros(n,3*n)];
L_6 = [zeros(n,varsigma+a+w+2*n),eye(n),zeros(n,2*n)];
L_7 = [zeros(n,varsigma+a+w+3*n),eye(n),zeros(n,n)];
L_8 = [zeros(n,varsigma+a+w+4*n),eye(n)];

Sigma_1 = [Pi_ux, Pi_um, zeros(m,w);
           zeros(w,n), zeros(w,a), eye(w)];
Lambda_1 = [L_8;L_2;L_3];
Lambda_2 = [eye(n),zeros(n,psi)];
mathcalG_1 = Lambda_2*A*L_1-Lambda_2*L_1+Lambda_2*B*Sigma_1*Lambda_1;
mathcalG_2 = C*L_1+D*Sigma_1*Lambda_1;
mathcalG_3 = [Pi_px*L_8+Pi_pm*L_2;L_2];

Upsilon_1d3 = Lambda_2*A*L_1+Lambda_2*B*Sigma_1*Lambda_1-L_4; % \Upsilon_{1,3}
Upsilon_2d3 = L_5-Lambda_2*A*L_1-Lambda_2*B*Sigma_1*Lambda_1; % \Upsilon_{2,3}
Upsilon_1 = [L_4;L_5;Upsilon_1d3;L_6+Upsilon_1d3];
Upsilon_2 = [-L_4;-L_5;Upsilon_2d3;L_7-L_5-Lambda_2*L_1];
Upsilon_3 = [L_4;L_5;zeros(n,bar_n);L_6];
Upsilon_4 = [L_4;L_5;zeros(n,bar_n);L_7];
Upsilon_5 = [zeros(n,bar_n);zeros(n,bar_n);Lambda_2*L_1-L_4;L_6-L_4];
Upsilon_6 = [zeros(n,bar_n);zeros(n,bar_n);L_5-Lambda_2*L_1;L_7-L_5];
Upsilon_7 = Upsilon_2-Upsilon_6;
Upsilon_8 = Upsilon_1-Upsilon_5;
Upsilon_9 = [Lambda_2*L_1-L_4;Lambda_2*L_1+L_4-L_6];
Upsilon_10 = [L_5-Lambda_2*L_1;L_5+Lambda_2*L_1-L_7];

epsilon_1 = 0.003; % This prescribed parameter embedded in Theorem can be adjusted for comparsion 0.3
epsilon_2 = 0.002; % This prescribed parameter embedded in Theorem can be adjusted for comparsion 0.2

vartheta_low = 1;   % The lower bound of sampling interval in ETC scheme
vartheta_upp = 5;  % The upper bound of sampling interval in ETC scheme 
                    % !!! Note that this value can be adjusted for comparsion
                    % while ensuring the feasibility of LMIs 
                    % baseline 设定值为5

vartheta_upp = 5; % 为了比较                  


cvx_begin sdp quiet
   cvx_solver mosek
   
   % Define variables
   variable P(varsigma,varsigma) symmetric
   variable T_1(n,n) symmetric
   variable T_2(n,n) symmetric
   variable Xi_1(n,n) symmetric
   variable Xi_2(n,n) symmetric
   variable R(4*n,4*n)
   variable N_1(bar_n,2*n)
   variable N_2(bar_n,2*n) 
   variable vgamma(a,1) % \gamma <--Variable name "gamma" is the name of a built-in MATLAB function.
                        %           Please choose a different name.
   P >= 1e-8*eye(varsigma); % P > 0
   T_1 >= 1e-8*eye(n); % T_1 > 0
   T_2 >= 1e-8*eye(n); % T_2 > 0
   Xi_1 >= 1e-8*eye(n); % Xi_1 > 0
   Xi_2 >= 1e-8*eye(n); % Xi_2 > 0
   
   M_gamma = diag(vgamma);
   mathcalM_DNN2 = [zeros(a,a), M_gamma; M_gamma zeros(a,a)];
   mathcalM_DNN = mathcalM_DNN1'*mathcalM_DNN2*mathcalM_DNN1; % \mathcal{M}_{DNN}
   
   mathcalQ = epsilon_1*L_4'*Xi_1*L_4 + epsilon_2*L_8'*Xi_1*L_8  -...
       (L_4-L_8)'*Xi_2*(L_4-L_8);
   amalg_1 = mathcalG_1'*T_2*mathcalG_1 + Upsilon_3'*R*Upsilon_7 +...
       Upsilon_7'*R'*Upsilon_3;
   amalg_2 = mathcalG_1'*T_1*mathcalG_1 + Upsilon_8'*R*Upsilon_4 +...
       Upsilon_4'*R'*Upsilon_8;
   
   mathcalG = (A*L_1+B*Sigma_1*Lambda_1)'*P*(A*L_1+B*Sigma_1*Lambda_1) -...
       mathcalG_1'*T_1*mathcalG_1 + mathcalG_1'*T_2*mathcalG_1 +...
       mathcalG_2'*M_Theta*mathcalG_2 + mathcalG_3'*mathcalM_DNN*mathcalG_3 +...
       (Upsilon_1'*R*Upsilon_2 - Upsilon_5'*R*Upsilon_6 + N_1*Upsilon_9 + N_2*Upsilon_10 ) +...
       (Upsilon_1'*R*Upsilon_2 - Upsilon_5'*R*Upsilon_6 + N_1*Upsilon_9 + N_2*Upsilon_10 )' -...
       L_1'*P*L_1 + mathcalQ;
   
   % Matrix inequalities
   [ mathcalG+vartheta_upp*amalg_1, vartheta_upp*N_1;...
       vartheta_upp*N_1', -vartheta_upp*blkdiag(T_1,3*T_1)] <= -1e-10*eye(bar_n+2*n);
   
   [ mathcalG+vartheta_upp*amalg_2, vartheta_upp*N_2;...
       vartheta_upp*N_2', -vartheta_upp*blkdiag(T_2,3*T_2)] <= -1e-10*eye(bar_n+2*n);
   
   for i = 1:a1
       [delta_rho^2 [W{1}(i,:),zeros(1,psi)];...
         [W{1}(i,:),zeros(1,psi)]', P] >= 0;  
   end
   
   [bar_varphi^2, [1,0,zeros(varsigma-n,1)];...
      [1,0,zeros(varsigma-n,1)]', P] >=0;
       
   minimize( trace(P(1:n,1:n)))
 
cvx_end
P; % To determine the inner approximation of robust RoA
Xi_1;
Xi_2;

if ~isnan(P)
    radii = 1./sqrt(eig(P(1:n,1:n)));
end

% figure (31)
% pvar x1 x2
% V = [x1,x2]*P(1:n,1:n)*[x1;x2];
% domain1 = [-10, 10, -10, 10];
% [C,h] = pcontour(V,1,domain1,'r');
% hold on

 %writematrix(P,'PlotROA/P1.csv') % 存储第一次实验数据 0.25 0.97
 %writematrix(P,'PlotROA/P2.csv') % 0.35 0.97
 %writematrix(P,'PlotROA/P3.csv') % 0.45 0.97

%writematrix(P,'PlotROA/P4.csv') % 0.25 1
% writematrix(P,'PlotROA/P5.csv') % 0.35 1  --->效果最好
%writematrix(P,'PlotROA/P6.csv') % 0.45 1

%writematrix(P,'PlotROA/P7.csv') % 0.25 1.03
%writematrix(P,'PlotROA/P8.csv') % 0.35 1.03
%writematrix(P,'PlotROA/P9.csv') % 0.45 1.03

%writematrix(P,'PlotROA/PP5.csv') % vartheta = 5;
%writematrix(P,'PlotROA/PP4.csv') % vartheta = 4;
%writematrix(P,'PlotROA/PP3.csv') % vartheta = 3;
%writematrix(P,'PlotROA/PP2.csv') % vartheta = 2;
%writematrix(P,'PlotROA/PP1.csv') % vartheta = 1;


Ellipsoid1=ellipsoid(inv(P(1:n,1:n)),zeros(n,1));
figure (333)
plot(Ellipsoid1,[1,2],'b');
hold on 
plot(0,0,'s');
hold on
hh = legend('$ \mathcal{A}=0.01 $', 'Equilibrium');
set(hh,'Interpreter','latex');



writematrix(P,sprintf('PlotROA/PP%d.csv', vartheta_upp)) 

%folderName = sprintf('PlotROA/P%d.csv', count);
%writematrix(P,folderName) 




