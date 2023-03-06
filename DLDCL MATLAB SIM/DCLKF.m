function [xi_EST] = DCLKF(N, G, rbt_ID, xj_Act, xj_IMU, xTg, mu_dx,IMU_cov,mu_cov)
% Kalman Fusion method for fusing multiple position estimates 
l_N = 1:N;
Neigh_ID = l_N(logical(G(rbt_ID,:)));

d_ji = zeros(length(xTg),N);

for j = 1:N
    d_ji(:,j) = RangeSensing(xj_Act(:,j),xj_Act(:,rbt_ID),mu_dx(j));
end
d_ji(:,rbt_ID) = zeros(length(xTg),1);

xihatj = zeros(length(xTg),length(Neigh_ID));
COVxij = zeros(length(xTg),length(xTg),length(Neigh_ID));

nj = 1;
for j = Neigh_ID
    xihatj(:,nj) = xj_IMU(:,j) + d_ji(:,j);
    COVxij(:,:,nj) = eye(length(xTg)).*(IMU_cov(j) + mu_cov(j)); 
    nj = nj + 1;
end

covM_sum = zeros(length(xTg),length(xTg));
fsum = zeros(length(xTg),1);
for k = 1:length(Neigh_ID)
      covM_sum = covM_sum + COVxij(:,:,k)^-1; 
      
      fsum = COVxij(:,:,k)\xihatj(:,k) + fsum;
end
xi_EST = covM_sum\fsum;

end

