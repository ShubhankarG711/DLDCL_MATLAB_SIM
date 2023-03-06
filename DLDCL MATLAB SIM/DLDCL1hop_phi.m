function [xi_Est,s_new,L_zi_new,Lpr_zi_new] = DLDCL1hop_phi(t, s_old, N, G, rbt_ID, xj_Act, xj_Prj, xj_IMU, xTg, mu_dx,L_wij,L_zi_prev,Lpr_wij,Lpr_zi_prev)
% DL-DCL algorithm for heading angle
To = 200; % periodic reset rate
eta_w = 2; % learning rate-1
eta_z = 2; % learning rate-2

% periodic reset
if floor(t/To) == ceil(t/To)
    L_wij = containers.Map();
    L_zi_prev = 0;
    
    Lpr_wij = containers.Map();
    Lpr_zi_prev = 0;
end

l_N = 1:N;
Neigh_ID = l_N(logical(G(rbt_ID,:))); % robot's neighbour set

% RPSS estimates 
for j = 1:N
    d_jH(:,j) = RelPhiSens(xj_Act(:,j),xTg,mu_dx(j));
    d_ji(:,j) = RelPhiSens(xj_Act(:,j),xj_Act(:,rbt_ID),mu_dx(j));
end
d_ji(:,rbt_ID) = zeros(length(xTg),1);

% get previous time cumulative loss values from containers 
nj = 1;
for j = Neigh_ID
   if isKey(L_wij,num2str(j))
       Lwij_prev(nj) = L_wij(num2str(j));
   else
       L_wij(num2str(j)) = 0;
       Lwij_prev(nj) = 0;
   end
   
   if isKey(Lpr_wij,num2str(j))
       Lwijpr_prev(nj) = Lpr_wij(num2str(j));
   else
       Lpr_wij(num2str(j)) = 0;
       Lwijpr_prev(nj) = 0;
   end
   
   nj = nj+1;
end

% calculate weights based on the previous time cumulative loss values
w_ij = exp(-eta_w.*Lwij_prev.*s_old)./sum(exp(-eta_w.*Lwij_prev.*s_old));
wp_ij = exp(-eta_w.*Lwijpr_prev.*s_old)./sum(exp(-eta_w.*Lwijpr_prev.*s_old));

zi = exp(-eta_z.*L_zi_prev.*s_old)./(exp(-eta_z.*L_zi_prev.*s_old) + exp(-eta_z.*Lpr_zi_prev.*s_old));

% form current time pose estimates of the landmark
nj = 1;
for j = Neigh_ID
   xtilde_jH(:,nj) = wrapToPi(xj_IMU(:,j)+d_jH(:,j));
   xhat_jH(:,nj) = wrapToPi(xj_Prj(:,j)+d_jH(:,j));
   nj = nj+1;
end
xtilde_H = wrapToPi(sum(repmat(w_ij,length(xTg),1).*xtilde_jH,2));
xhat_H = wrapToPi(sum(repmat(wp_ij,length(xTg),1).*xhat_jH,2));

xH_Est = wrapToPi(zi*xtilde_H + (1-zi)*xhat_H);

% loss function definition
lssfn = @(x,y) min(abs(wrapToPi(x-y))./(15*pi/180),1);

% update the cumulative loss values for the current time step
Lwij = (Lwij_prev.*s_old + lssfn(xtilde_jH,xTg))./(s_old + 1);
Lwijpr = (Lwijpr_prev.*s_old + lssfn(xhat_jH,xTg))./(s_old + 1);

L_zi = (L_zi_prev.*s_old + lssfn(xtilde_H,xTg))./(s_old + 1);
Lpr_zi = (Lpr_zi_prev.*s_old + lssfn(xhat_H,xTg))./(s_old + 1);

% update weights based on the current time cumulative loss values 
w_ij = exp(-eta_w.*Lwij.*(s_old+1))./sum(exp(-eta_w.*Lwij.*(s_old+1)));
wp_ij = exp(-eta_w.*Lwijpr.*(s_old+1))./sum(exp(-eta_w.*Lwijpr.*(s_old+1)));

zi = exp(-eta_z.*L_zi.*(s_old+1))./(exp(-eta_z.*L_zi.*(s_old+1)) + exp(-eta_z.*Lpr_zi.*(s_old+1)));

% change H to i : (Landmark to robot i)
% form current time pose estimates of i^th robot (self)
nj = 1;
for j = Neigh_ID
   xtilde_ji(:,nj) = wrapToPi(xj_IMU(:,j)+d_ji(:,j));
   xhat_ji(:,nj) = wrapToPi(xj_Prj(:,j)+d_ji(:,j));
   nj = nj+1;
end
xtilde_i = wrapToPi(sum(repmat(w_ij,length(xTg),1).*xtilde_ji,2));
xhat_i = wrapToPi(sum(repmat(wp_ij,length(xTg),1).*xhat_ji,2));

xi_Est = wrapToPi(zi*xtilde_i + (1-zi)*xhat_i);

% Update Cumulative loss value containers
nj = 1;
for j = Neigh_ID
   L_wij(num2str(j)) = Lwij(nj);
   Lpr_wij(num2str(j)) = Lwijpr(nj);
   
   nj = nj + 1;
end
L_zi_new = L_zi;
Lpr_zi_new = Lpr_zi;

% update counter 
if floor(t/To) == ceil(t/To)
    s_new = 0;  
else
    s_new = s_old + 1;
end
