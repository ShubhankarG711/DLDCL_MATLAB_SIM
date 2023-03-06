function [x_next,phi_next,delta_x,delta_phir] = robotctrldyNnoise(dT_safe, kvTf_i, kvNs_i, kvRT_i, kw_i, deltaT, x_curr, phi_curr, RM_curr, x_nrst_prev, x_nrst_nxt, x_nrst_rtb, x_tg_curr, x_noise, ID_prev, ID_nxt, ID_rtb, x_curr_act, phi_curr_act)
% robot kinematics with its control law

phi_curr = wrapToPi(phi_curr);
phi_curr_act = wrapToPi(phi_curr_act);

delta_x_T = x_tg_curr - x_curr;

delta_x_N_prev = x_nrst_prev - x_curr;
delta_x_N_nxt = x_nrst_nxt - x_curr;
delta_x_N_rtb = x_nrst_rtb - x_curr;

vbar_i = kvTf_i*RM_curr'*(delta_x_T/norm(delta_x_T))*(norm(delta_x_T)-dT_safe);

delta_v_N = -kvNs_i*RM_curr'*(delta_x_N_rtb/(norm(delta_x_N_rtb))^2);   
delta_v_N = min(max(delta_v_N,-5),5);

vbar_rot3d = kvRT_i*cross([0 0 1]',[-delta_x_T' 0]');
vbar_rot = RM_curr'*vbar_rot3d(1:2,1);

headTar = atan2(-delta_x_T(2),-delta_x_T(1));

headTg_eff = 1.0*headTar;

delta_phi_target = wrapToPi(headTg_eff - phi_curr);

v_ctrl = min(max(vbar_i + delta_v_N + vbar_rot,-10),10);
w_ctrl = min(max(kw_i*delta_phi_target,-45*pi/180),45*pi/180);

x_next = x_curr_act + deltaT*RM_curr*v_ctrl + x_noise.*[1;1].*0; % zero disturbance
phi_next = wrapToPi(phi_curr_act + deltaT*w_ctrl);

delta_x = deltaT*RM_curr*v_ctrl;
delta_phir = deltaT*w_ctrl;
end