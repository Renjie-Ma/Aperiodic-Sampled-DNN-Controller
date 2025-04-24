%%%%%%%%% ETDNN_Stability
%%%%%%%%% Renjie Ma, Harbin Institute of Technology
%%%%%%%%% Dec 2023
clear;
clc;

%% System Parameters
Ts = 0.01; 
m = 0.15;
l = 0.5;
mu = 0.05;
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
x0 =  [0.23; 3.5];
 % 第一次测试：[0.23; 3.5] ； 第二次测试： [0.43 ;-3.0] ； 第三次测试：
%[-0.33; -3.3];
%x0 = [0.43;-3.0];
eta0 = [x0;xi0];
Nstep = 800;

% mented state updates of Time-triggered DNN-controlled system 
[eta,u] = T2NNClosedLoop(Nstep,eta0,W,bias,A,B,umax,umin,statenum); % Adjust the saturation via T2NNClosedLoop function
% save the state/control data pairs
writematrix(eta,'t2eta.csv')  
writematrix(u,'t2u.csv')
    
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

delta_rho = 0.4; % This value can be adjusted for comparsion 
delta_sigma = 1; % This value can be adjusted for comparsion
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


Ellipsoid1=ellipsoid(inv(P(1:n,1:n)),zeros(n,1));
figure (3)
plot(Ellipsoid1,[1,2],'b');
hold on 
plot(0,0,'s');
hold on
hh = legend('$ \mathcal{A}=0.01 $', 'Equilibrium');
set(hh,'Interpreter','latex');

%% ETC Scheme

vartheta =1; % This value belongs to [1,vartheta_upp]

mu_et = 0.2; % \mu =0.02 0.2
gg =480;  % g=2 480
precon = 1-mu_et-(1/gg); % 1-\mu-g^{-1}>0

if precon >0 
    label_etc = 1; % The prerequisite for ETC is satisfied
else
    label_etc = 0; % The parameters of \mu and g need to be re-chosen
end

T = vartheta; % sampling interval: vartheta_low <=vartheta<= vartheta_upp; This parameter can be adjusted for comparsion
pts = 1:T:Nstep; % sampling instants
nS = length(pts);    % calculate the number of sampling instants

isi = 0;
iE = 1;
seqi(1) = 1;

alpha(1) = 5000000; % The initial value of dynamic variable embedded in ETC
cxS{1} = x0; % The initial value of sampler
cxE{1} = x0; % The initial value of event trigger 

cx{1} = x0; % Revisit the intial state x(0)
etaa{1} = eta0; % Revisit the intial state \eta(0)

for k = 1: Nstep
    time(k) = k-1;
    
    for s= 1 : nS  % sampling times nS  = Nstep/vartheta 
        si = pts(s); % sampling instants si 
            
        %tjq: \theta{j,q} = k_{q} + j*vartheta (vartheta_upp)
        % \theta_{j,q}: 当前采样点，前后两个采样点间距长度为 vartheta
        % 对于一段采样间隔，长度为vartheta
        if k == si  % si 指代采样点！！事件触发点，初始时刻假定触发，循环更新
           isi = isi +1; % counting the number of sampling intervals, each with the length vartheta_upp 
                         % 开启第一段采样间隔[0,vartheta]
          
           cxS{isi} = cx{si}; % cxS{isi} = x(\theta_{j,q}) the current-sampled state 
                              % 将当前采样状态赋给cxS{isi}（对应第isi次采样，采样过程中间状态保持）; 
           beta(isi) = epsilon_1*(cxS{isi}-x_eq)'*Xi_1*(cxS{isi}-x_eq) +...
                       epsilon_2*(cxE{iE}-x_eq)'*Xi_1*(cxE{iE}-x_eq) -... 
                       (cxS{isi}-cxE{iE})'*Xi_2*(cxS{isi}-cxE{iE});  % cxE{iE} = x(k_{q}); iE 为触发时刻，记作k_q
           alpha(isi+1) = (1-mu_et)*alpha(isi) + beta(isi);  % ！！ alpha(ist+1) 对应着下一个采样间隔,间距vartheta  
                  
        
            % 判断采样时刻是否为事件触发时刻
            if alpha(isi)+ gg* beta(isi) <0
                   %  满足事件触发条件
               iE =iE + 1; % counting the number of triggering instants 计数事件触发次数
               seqi(iE) = si; % 抽取当前采样点构成采样序列
               intiE(iE) = seqi(iE)-seqi(iE-1); % 计算同前一次采样的时间间隔
               cxE{iE} = cxS{isi}; % 对触发状态进行更新
            end 
            
        end
                     
    end
    
    x = cxE{iE};
    uuu = NN_Controller(x,W,bias,umax,umin);
    
    uuetc{iE} = uuu;
    
    
       etaa{k+1} = A*etaa{k}+B*[uuetc{iE}; etaa{k}(1,:)-sin(etaa{k}(1,:))];
       cx{k+1} = etaa{k+1}(1:2,:);
    
    
    etaa1(k) = etaa{k}(1);  
    etaa2(k) = etaa{k}(2);
    etaa3(k) = etaa{k}(3);
end

folderName = sprintf('%.2f_%.2f', x0(1), x0(2));
folderName = strrep(folderName, '-', 'm');
if ~isfolder(folderName)
    mkdir(folderName); % 如果不存在，则创建   
end



if vartheta==4
    %第一次测试数据存储： vartheta = 4
    writematrix(etaa1,fullfile(folderName, 'etaa11.csv'))  % 1etaa --->checkepsilon_1 = 0.8
    writematrix(etaa2,fullfile(folderName, 'etaa12.csv'))            % checkepsilon_2 = 0.6
    writematrix(etaa3,fullfile(folderName, 'etaa13.csv'))            % bar_s = 10
    % writematrix(time,'time')                   % 0.23;3.5
    writematrix(seqi,fullfile(folderName, 'seqi1.csv'))
    writematrix(intiE,fullfile(folderName, 'intiE1.csv'))
end

if vartheta==3
    %第二次测试数据存储： vartheta = 3
    writematrix(etaa1,fullfile(folderName, 'etaa21.csv'))  % 1etaa --->checkepsilon_1 = 0.8
    writematrix(etaa2,fullfile(folderName, 'etaa22.csv'))            % checkepsilon_2 = 0.6
    writematrix(etaa3,fullfile(folderName, 'etaa23.csv'))            % bar_s = 10
    %writematrix(time,'time')                   % 0.23;3.5
    writematrix(seqi,fullfile(folderName, 'seqi2.csv'))
    writematrix(intiE,fullfile(folderName, 'intiE2.csv'))
end

if vartheta==2
%第三次测试数据存储： vartheta = 2
  writematrix(etaa1,fullfile(folderName, 'etaa31.csv'))  % 1etaa --->checkepsilon_1 = 0.8
  writematrix(etaa2,fullfile(folderName, 'etaa32.csv'))            % checkepsilon_2 = 0.6
  writematrix(etaa3,fullfile(folderName, 'etaa33.csv'))            % bar_s = 10
  %writematrix(time,'time')                   % 0.23;3.5
  writematrix(seqi,fullfile(folderName, 'seqi3.csv'))
  writematrix(intiE,fullfile(folderName, 'intiE3.csv'))
end

if vartheta==1
    %第四次测试数据存储： vartheta = 1
   writematrix(etaa1,fullfile(folderName, 'etaa41.csv'))  % 1etaa --->checkepsilon_1 = 0.8
   writematrix(etaa2,fullfile(folderName, 'etaa42.csv'))            % checkepsilon_2 = 0.6
   writematrix(etaa3,fullfile(folderName, 'etaa43.csv'))            % bar_s = 10
    %writematrix(time,'time')                   % 0.23;3.5
   writematrix(seqi,fullfile(folderName, 'seqi4.csv'))
   writematrix(intiE,fullfile(folderName, 'intiE4.csv'))
end

figure(44)
subplot(3,1,1)
plot(Ts*time,etaa1,'-','linewidth',1.0);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('$x_{1}$');
set(hh, 'Interpreter', 'latex');
subplot(3,1,2)
plot(Ts*time,etaa2,'-','linewidth',1.0);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('$x_{1}$');
set(hh, 'Interpreter', 'latex');
subplot(3,1,3)
plot(Ts*time,etaa3,'-','linewidth',1.0);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('$x_{1}$');
set(hh, 'Interpreter', 'latex');

figure(114)

stem(seqi,intiE,'color','(0.0,0.45,0.74)','linewidth',1.0);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex','fontsize',13); 
% yl=ylabel('Inter-event intervals');
yl=ylabel('$d_{j}$');
set(yl,'Interpreter','latex','fontsize',13); 
title('\fontsize{13}(a)') 
%axis([15,200,0,11])
% set(gca,'FontSize',12) 


















































