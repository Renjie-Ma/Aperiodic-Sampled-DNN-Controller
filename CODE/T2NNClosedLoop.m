function [eta,u] = T2NNClosedLoop(Nstep,eta0,W,bias,A,B,umax,umin,statenum)


% This function aims at describing the time-triggered augmented dynamics Eq.(8)

u0 = NN_Controller(eta0(1:statenum,:),W,bias,umax,umin);
omega0 = eta0(1,:)-sin(eta0(1,:));
baru0 = [u0,omega0];
Nbaru = numel(baru0); % return the number of elements embedded in baru0
Nu = numel(u0); % return the number of elements embedded in u0
Nx = numel(eta0(1:statenum,:)); % return the number of elements embedded in x0
Neta = numel(eta0); % return the number of elements embedded in eta0
Nomega = numel(omega0); % return the number of elements embedded in omega0

eta = zeros(Neta,Nstep);
eta(:,1) = eta0;

u = zeros(Nu,Nstep);
u(:,1) = u0;

baru = zeros(Nbaru,Nstep);
baru(:,1) = baru0;

omega = zeros(Nomega,Nstep);
omega(:,1) = omega0;

% Simulate the augmented plant

for i = 2: Nstep
    
    eta(:,i) = A*eta(:,i-1)+B*baru(:,i-1);
    u(:,i) = NN_Controller(eta(1:2,i),W,bias,umax,umin);
    omega(:,i) = eta(1,i)-sin(eta(1,i));
    baru(:,i) = [u(:,i);omega(:,i)];
    
end

end




