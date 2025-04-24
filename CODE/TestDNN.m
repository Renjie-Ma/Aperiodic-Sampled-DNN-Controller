function [x,u] = TestDNN(Nstep,x0,W,bias,umax,umin,statenum)


% This function aims at describing the time-triggered augmented dynamics Eq.(8)

u0 = NN_Controller(x0,W,bias,umax,umin);
Nu = numel(u0); % return the number of elements embedded in u0
Nx = numel(x0); % return the number of elements embedded in x0


u = zeros(Nu,Nstep);
u(:,1) = u0;

x = zeros(Nx,Nstep);
x(:,1) = x0;

Ts = 0.02; 
m = 0.15;
l = 0.5;
mu = 0.05;
g = 9.8;

A = [1 Ts; 0 1-Ts*mu/(m*l*l)];
B = [0; Ts/(m*l*l)];
F = [0; Ts*g/l];


% Simulate the augmented plant

for i = 2: Nstep
    
    x(:,i) = A*x(:,i-1)+B*u(:,i-1)+F*sin(x(1,i-1));
    u(:,i) = NN_Controller(x(:,i),W,bias,umax,umin);
    
    
end

end
