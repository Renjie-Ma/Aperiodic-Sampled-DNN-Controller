function [A,B,C,D,statenum,colnum_A_Phi,omeganum,M_Theta] = Pendulum(Ts,m,l,mu,g,ls,ms)

% System parameters of inverted pendulum with IQC virtual filter

A_Gamma = [1, Ts; 
           Ts*g/l, 1-Ts*mu/(m*l*l)];
B_Gamma = [0; 
           Ts/(m*l*l)];
F_Gamma = [0; 
           -Ts*g/l];
C_Gamma = [1, 0];
D_Gamma = 0; 
G_Gamma = 0; 

rownum_A_Gamma = size(A_Gamma,1); % return the row number of matrix A_Gamma
statenum = rownum_A_Gamma;
omeganum = size(F_Gamma,2);

% Define the virtual filter by off-by-one IQC [Lessard2016Analysis,SIAM]
A_Phi = 0; 
B_Phi = -ls;
F_Phi = 1;
C_Phi = [1; 
         0];
D_Phi = [ls;
         -ms];
G_Phi = [-1;
         1];
colnum_A_Phi = size(A_Phi,2); % return the column number of matrix A_Phi
M_Theta = [0,1;
           1,0];

A = [A_Gamma, zeros(rownum_A_Gamma,colnum_A_Phi);
    B_Phi*C_Gamma, A_Phi];
B = [B_Gamma, F_Gamma;
    B_Phi*D_Gamma, B_Phi*G_Gamma+F_Phi];
C = [D_Phi*C_Gamma,C_Phi];
D = [D_Phi*D_Gamma, D_Phi*G_Gamma+G_Phi];



end

