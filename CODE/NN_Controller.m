function u = NN_Controller(x,W,bias,umax,umin)

% This function aims at formulating the DNN controller

layernumber = numel(W);
pp = x;
for i = 1: (layernumber-1)
    pp = W{i}*pp+bias{i};
    pp = tanh(pp);
end

u = W{end}*pp+bias{end};
u = SAT(u,umax,umin); % The saturation 

end


