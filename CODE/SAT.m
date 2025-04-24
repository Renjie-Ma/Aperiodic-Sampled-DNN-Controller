function usat = SAT(u,umax,umin)

% This function aims at formulation the saturation of control input

if u > umax
    usat = umax;
elseif umin <=u && u <=umax
    usat = u;
else
    usat = umin;
end
