function [x_next, phi_next] = targetctrldyn(deltaT, x_curr, phi_curr, RM_curr, vbar_curr, wbar_curr)
% 3-DOF kinematic model for the target/Landmark pose simulation 

phi_curr = wrapToPi(phi_curr);

x_next = x_curr + deltaT*RM_curr*vbar_curr;
phi_next = wrapToPi(phi_curr + deltaT*wbar_curr);

end

