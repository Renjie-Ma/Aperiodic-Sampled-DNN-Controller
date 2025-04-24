%%%%%%%%% STDNN_Stability
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
umax =0.9; 
umin = -umax; % saturation

%% Define the initial values

% augmented state eta(0)=eta0
xi0 = zeros(colnum_A_Phi,1);
%x0 = [0.23;3.5]; % case 1

%x0 = [0.43; 3.0]; % case 2

%x0 = [-0.33;-1.8]; %case 3 change into [-0.33;-3.3]
x0 = [0.19;3.5];


eta0 = [x0;xi0];
Nstep = 800;
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
tilde_n = varsigma+a+w+4*n;
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

delta_rho = 0.5; % This value can be adjusted for comparsion 
delta_sigma = 1.025; % This value can be adjusted for comparsion

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

%% Linear Matrix Inequality (c.t.: Svent-Triggered Case(Case 3))

overlineL_1 = [eye(varsigma),zeros(varsigma,tilde_n-varsigma)]; 
overlineL_2 = [zeros(a,varsigma), eye(a), zeros(a,tilde_n-varsigma-a)];
overlineL_3 = [zeros(w,varsigma+a), eye(w),zeros(w,4*n)];
overlineL_4 = [zeros(n,varsigma+a+w), eye(n), zeros(n,3*n)];
overlineL_5 = [zeros(n,varsigma+a+w+n), eye(n), zeros(n,2*n)];
overlineL_6 = [zeros(n,varsigma+a+w+2*n), eye(n), zeros(n,n)];
overlineL_7 = [zeros(n,varsigma+a+w+3*n), eye(n)];

Sigma_1 = [Pi_ux, Pi_um, zeros(m,w);
           zeros(w,n), zeros(w,a), eye(w)];
       
overlineLambda_1 = [overlineL_4;overlineL_2;overlineL_3];
Lambda_2 = [eye(n),zeros(n,psi)];
overlinemathcalG_1 = Lambda_2*A*overlineL_1 - Lambda_2*overlineL_1 + Lambda_2*B*Sigma_1*overlineLambda_1;
overlinemathcalG_2 = C*overlineL_1 + D*Sigma_1*overlineLambda_1;
overlinemathcalG_3 = [Pi_px*overlineL_4+Pi_pm*overlineL_2; overlineL_2];

overlineUpsilon_12 = Lambda_2*A*overlineL_1 + Lambda_2*B*Sigma_1*overlineLambda_1; % \overline{Upsilon}_{12}
overlineUpsilon_13 = overlineL_5 - Lambda_2*A*overlineL_1 - Lambda_2*B*Sigma_1*overlineLambda_1; % \overline{Upsilon}_{13}

overlineUpsilon_1 = [ overlineL_4; overlineL_5; overlineUpsilon_12; overlineL_6+overlineUpsilon_12]; % \overline{Upsilon}_{1}
overlineUpsilon_2 = [-overlineL_4; -overlineL_5; overlineUpsilon_13; overlineL_7-overlineL_5-Lambda_2*overlineL_1]; % \overline{Upsilon}_{2}
overlineUpsilon_3 = [overlineL_4; overlineL_5; zeros(n,tilde_n); overlineL_6];  % \overline{Upsilon}_{3} 
overlineUpsilon_4 = [overlineL_4; overlineL_5; zeros(n,tilde_n); overlineL_7]; % \overline{Upsilon}_{4}
overlineUpsilon_5 = [zeros(n,tilde_n); zeros(n,tilde_n); Lambda_2*overlineL_1-overlineL_4; overlineL_6-overlineL_4]; % \overline{Upsilon}_{5}
overlineUpsilon_6 = [zeros(n,tilde_n); zeros(n,tilde_n); overlineL_5-Lambda_2*overlineL_1; overlineL_7-overlineL_5]; % \overline{Upsilon}_{6}
overlineUpsilon_7 = overlineUpsilon_2-overlineUpsilon_6; % \overline{Upsilon}_{7}
overlineUpsilon_8 = overlineUpsilon_1-overlineUpsilon_5; % \overline{Upsilon}_{8}
overlineUpsilon_9 = [Lambda_2*overlineL_1-overlineL_4; Lambda_2*overlineL_1+overlineL_4-overlineL_6]; % \overline{Upsilon}_{9}
overlineUpsilon_10 = [overlineL_5-Lambda_2*overlineL_1; overlineL_5+Lambda_2*overlineL_1-overlineL_7]; % \overline{Upsilon}_{10}
overlineUpsilon_11 = [overlineL_5; overlineL_4]; % \overline{Upsilon}_{11}

checkepsilon_1 = 0.8; % \check{epsilon}_{1} embedded in Theorem can be adjusted for comparsion 0.8
checkepsilon_2 = 0.6; % \check{epsilon}_{2} embedded in Theorem can be adjusted for comparsion 0.6
bar_s = 8; % The upper bound of triggering interval in STC scheme > 1 100

cvx_begin sdp quiet
   cvx_solver mosek
   
   % Define variables
   variable overlineP(varsigma,varsigma) symmetric   % \overline{P}
   variable overlineT_1(n,n) symmetric               % \overline{T}_{1}
   variable overlineT_2(n,n) symmetric               % \overline{T}_{2}
   variable checkXi_1(n,n) symmetric                 % \check{Xi}_{1}
   variable checkXi_2(n,n) symmetric                 % \check{Xi}_{2}
   variable overlineR(4*n,4*n)                       % \overline{R}
   variable overlineN_1(tilde_n,2*n)                 % \overline{N}_{1}
   variable overlineN_2(tilde_n,2*n)                 % \overline{N}_{2}
   variable vgamma(a,1)                              % \gamma_{i}>=0
   variable barlambda_1(1,1)                         % \bar{lambda}_{1} 
   variable barlambda_2(1,1)                         % \bar{lambda}_{2}
   
   variable MM_Theta(n,n) symmetric                   % if treating M_Theta as a decision variable
                        
   overlineP >= 1e-8*eye(varsigma); % \overline{P} > 0
   overlineT_1 >= 1e-10*eye(n); % \overline{T}_{1} > 0
   overlineT_2 >= 1e-10*eye(n); % \overline{T}_{2}
   checkXi_1 >= 1e-8*eye(n); % \check{Xi}_{1} > 0
   checkXi_2 >= 1e-8*eye(n); % \check{Xi}_{2} > 0
   MM_Theta >= 1e-8*eye(n);
   
   M_gamma = diag(vgamma);
   mathcalM_DNN2 = [zeros(a,a), M_gamma; M_gamma zeros(a,a)];
   mathcalM_DNN = mathcalM_DNN1'*mathcalM_DNN2*mathcalM_DNN1; % \mathcal{M}_{DNN}
   
   % Define the matrix \tilde{\mathcal{G}}
   
   mathcalX_1 = [checkepsilon_1*checkXi_1-checkXi_2; checkXi_2']; % \mathcal{X}_{1}
   mathcalX_2 = [checkXi_2; checkepsilon_2*checkXi_1-checkXi_2];  % \mathcal{X}_{2}
   widetildeXi = [mathcalX_1,mathcalX_2];
     
   tildemathcalG = (A*overlineL_1+B*Sigma_1*overlineLambda_1)'*overlineP*(A*overlineL_1+B*Sigma_1*overlineLambda_1) +...
       (bar_s -1)*overlinemathcalG_1'*overlineT_1*overlinemathcalG_1 +...
       (bar_s +1)*overlinemathcalG_1'*overlineT_2*overlinemathcalG_1 +...
       bar_s*barlambda_1*overlineUpsilon_7'*overlineUpsilon_7 + bar_s*barlambda_2*overlineUpsilon_4'*overlineUpsilon_4 +...
       overlinemathcalG_2'*MM_Theta*overlinemathcalG_2 + overlinemathcalG_3'*mathcalM_DNN*overlinemathcalG_3 +...
       (overlineUpsilon_1'*overlineR*overlineUpsilon_2 - overlineUpsilon_5'*overlineR*overlineUpsilon_6 +...
       overlineN_1*overlineUpsilon_9 + overlineN_2*overlineUpsilon_10) +...
       (overlineUpsilon_1'*overlineR*overlineUpsilon_2 - overlineUpsilon_5'*overlineR*overlineUpsilon_6 +...
       overlineN_1*overlineUpsilon_9 + overlineN_2*overlineUpsilon_10)' -...
       overlineUpsilon_11'*widetildeXi*overlineUpsilon_11 - overlineL_1'*overlineP*overlineL_1; 
       % MM_Theta: if treat M_Theta as a decision variable
       
   % Define the matrix \eth_{1}
   eth_1 = [bar_s*overlineN_1, bar_s*overlineN_2, bar_s*overlineUpsilon_3'*overlineR, bar_s*overlineUpsilon_8'*overlineR];
   
   % Define the matrix \eth_{2}
   overlinemathcalT_1 = blkdiag(overlineT_1, 3*overlineT_1);
   overlinemathcalT_2 = blkdiag(overlineT_2, 3*overlineT_2);
   
   eth_2 = blkdiag(bar_s*overlinemathcalT_1,bar_s*overlinemathcalT_2,bar_s*barlambda_1*eye(4*n), bar_s*barlambda_2*eye(4*n) );
   
   % Matrix inequalities
   [ tildemathcalG, eth_1;...
       eth_1', eth_2] <= -1e-10*eye(tilde_n+2*n+2*n+4*n+4*n);
   
   for i = 1:a1
       [delta_rho^2 [W{1}(i,:),zeros(1,psi)];...
         [W{1}(i,:),zeros(1,psi)]', overlineP] >= 0;  
   end
   
   [bar_varphi^2, [1,0,zeros(varsigma-n,1)];...
      [1,0,zeros(varsigma-n,1)]', overlineP] >=0;
  
  minimize( trace(overlineP(1:n,1:n)))
   
cvx_end
overlineP;
checkXi_1;
checkXi_2;

% Ellipsoid1=ellipsoid(inv(overlineP(1:n,1:n)),zeros(n,1)); 
% figure (114)
% plot(Ellipsoid1,[1,2],'b');
% hold on 
% plot(0,0,'s');
% hold off
% hh = legend('$ \mathcal{A}=0.01 $', 'Equilibrium');
% set(hh,'Interpreter','latex');

%% STC scheme

isi = 1; % 自触发次数计数

seqi(1) = 0; % 初始化自触发时刻序列
intiE(1) =0; % 初始化自触发间隔

cxE{1} = x0; % The initial value of event trigger 

cx{1} = x0; % Revisit the intial state x(0)
etaa{1} = eta0; % Revisit the intial state \eta(0)

si =1; % 初始自触发时刻


for k = 1: Nstep % 时域上界为Nstep
    time(k) = k-1; 
    
     if k == si  % si 指代自触发时刻； 结尾需要更新
        etaaaS{isi} = etaa{si}; % 提取自触发时刻的系统状态; 结尾需要更新 etaa{si}---> eta(k_{q})
         
        count_q = 0; % 引入for 循环计数变量
        etaaa{1} = etaaaS{isi}; % 用以代入for循环计算 \eta{k_q+s_k} 
        
        for s = 1: bar_s % bar_s 为设定的最大自触发时刻间隔 -->该循环用以决定s_k的值，该循环中，实际上k ==si， si 为上一个自触发时刻
           
            x= etaaaS{isi}(1:2,:); % 设置控制输入为最近一次自触发时刻的状态
            uu = NN_Controller(x,W,bias,umax,umin); 
            
            etaaa{s+1} = A*etaaa{s}+B*[uu; etaaa{s}(1,:)-sin(etaaa{s}(1,:)) ]; % 用于计算 eta(k_q+s_k)
            
            check_ST(s) = [etaaa{s+1}(1:2,:)-x_eq; etaaaS{isi}(1:2,:)-x_eq]' *...
                       [checkepsilon_1*checkXi_1-checkXi_2, checkXi_2;
                        checkXi_2',checkepsilon_2*checkXi_1-checkXi_2] *...
                       [etaaa{s+1}(1:2,:)-x_eq; etaaaS{isi}(1:2,:)-x_eq];  % 计算 \mathcal{S}(x(k_q),s_k)
            count_q = count_q+1; % 计数 for 循环迭代次数
                         
            %判断当前时刻是否为自触发时刻
            if check_ST(s) >= 0
                %|| count_q == bar_s
                ss = count_q;  % 下一个自触发时刻为 si+ss
                % break % 注释掉break 是因为想要找到尽可能最大的ss, 如果有break的话，最小值就会跳出
            end
            
        end
        
        si = si + ss; % 更新下一个自触发时刻
        isi =isi +1; % 自触发次数加1
       
        seqi(isi) = si; % 抽取当前时刻构成自触发序列
        intiE(isi) = seqi(isi)-seqi(isi-1); % 计算同前一次自触发时刻的事件间隔

       % cxE{isi} = etaa{si}(1:2,:); % 更新自触发中的x(k_{q+1})        
     end
     
     
     if seqi(isi-1) <= k <seqi(isi)
         
       x= etaa{si-ss}(1:2,:); % 设置控制输入为最近一次自触发时刻的状态
       uu = NN_Controller(x,W,bias,umax,umin); 
       etaa{k+1} = A*etaa{k}+B*[uu; etaa{k}(1,:)-sin(etaa{k}(1,:))];
       cx{k+1} = etaa{k+1}(1:2,:);
       
     end

    
    etaa1(k) = etaa{k}(1);  
    etaa2(k) = etaa{k}(2);
    etaa3(k) = etaa{k}(3);
    
end

% 第一次测试数据存储： 
% writematrix(etaa1,'etaa11.csv')  % 1etaa --->checkepsilon_1 = 0.8
% writematrix(etaa2,'etaa12.csv')            % checkepsilon_2 = 0.6
% writematrix(etaa3,'etaa13.csv')            % bar_s = 10
% writematrix(time,'time')                   % 0.23;3.5
% writematrix(seqi,'seqi1.csv')
% writematrix(intiE,'1intiE1.csv')

% 第二次测试数据存储： 
% writematrix(etaa1,'etaa21.csv')  % 1etaa --->checkepsilon_1 = 0.8
% writematrix(etaa2,'etaa22.csv')             % 0.43;-3.0
% writematrix(etaa3,'etaa23.csv')
% %writematrix(time,'time')
% writematrix(seqi,'seqi2.csv')
% writematrix(intiE,'intiE2.csv')

folderName = sprintf('%.2f_%.2f', x0(1), x0(2));
folderName = strrep(folderName, '-', 'm');
if ~isfolder(folderName)
    mkdir(folderName); % 如果不存在，则创建
end

% 第三次测试数据存储： 
writematrix(etaa1,fullfile(folderName, 'staa1.csv')) % 1etaa --->
writematrix(etaa2,fullfile(folderName, 'staa2.csv'))
writematrix(etaa3,fullfile(folderName, 'staa3.csv'))
%writematrix(time,'time')
writematrix(seqi,fullfile(folderName, 'STseqi.csv'))
writematrix(intiE,fullfile(folderName, 'STintiE.csv'))

figure(55)
subplot(3,1,1)
plot(Ts*time,etaa1,'-','linewidth',1.5);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('Self-Triggered Neural-Feedbacked $x_{1}(k)$','FontSize',10);
set(hh, 'Interpreter', 'latex');
title('\fontsize{11}(a) The state trajectory of angular position') 
subplot(3,1,2)
plot(Ts*time,etaa2,'-','linewidth',1.5);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('Self-Triggered Neural-Feedbacked $x_{2}(k)$','FontSize',10);
set(hh, 'Interpreter', 'latex');
title('\fontsize{11}(b) The state trajectory of angular velocity') 
subplot(3,1,3)
plot(Ts*time,etaa3,'-','linewidth',1.5);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex');
hh = legend('Self-Triggered Neural-Feedbacked $\xi(k)$','FontSize',10);
set(hh, 'Interpreter', 'latex');

figure(115)

stem(seqi,intiE,'color','(0.0,0.45,0.74)','linewidth',1.0);
xl=xlabel('Time (Sec.)');
set(xl,'Interpreter','latex','fontsize',13); 
% yl=ylabel('Inter-event intervals');
yl=ylabel('$d_{j}$');
set(yl,'Interpreter','latex','fontsize',13); 
%axis([400,600,0,11])
% title('\fontsize{13}(a)') 
% set(gca,'FontSize',12) 




















