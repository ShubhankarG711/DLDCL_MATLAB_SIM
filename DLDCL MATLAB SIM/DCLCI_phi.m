function [xi_EST] = DCLCI_phi(N, G, rbt_ID, xj_Act, xj_IMU, xTg, mu_dx,IMU_cov,mu_cov)
% Covariance Intersection method for fusing multiple yaw angle estimates
l_N = 1:N;
Neigh_ID = l_N(logical(G(rbt_ID,:)));

d_ji = zeros(length(xTg),N);
for j = 1:N
    d_ji(:,j) = RelPhiSens(xj_Act(:,j),xj_Act(:,rbt_ID),mu_dx(j));
end
d_ji(:,rbt_ID) = zeros(length(xTg),1);

xihatj = zeros(length(xTg),length(Neigh_ID));
COVxij = zeros(length(xTg),length(xTg),length(Neigh_ID));

nj = 1;
for j = Neigh_ID
    xihatj(:,nj) = wrapToPi(xj_IMU(:,j) + d_ji(:,j));
    COVxij(:,:,nj) = eye(length(xTg)).*(IMU_cov(j) + mu_cov(j)); 
    nj = nj + 1;
end

if length(Neigh_ID) > 1
      [xi_EST,~] = fusecovint(xihatj,COVxij);
       xi_EST = wrapToPi(xi_EST);
else
       xi_EST = xihatj;
end

end