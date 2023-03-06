%%% DL-DCL Matlab Scalability Study Master Code; comparative study with KF, CI, CU %%%
%%% by Shubhankar Gupta, Suresh Sundaram %%%
%%% AIRL, Dept. of Aerospace Engg., IISc. Bangalore, India %%%
%%% For more details, please check the associated AAAI-23 conference paper titled,
%%% 'Moving-Landmark assisted Distributed Learning based Decentralized
%%% Cooperative Localization (DL-DCL) with Fault Tolerance.' %%%

%%% This code double loops over a modified version of the main code in the file 
%%% for Matlab Simulation Video and generates scalability study plots with 
%%% results averaged over multiple simulation runs. 
%%% For detailed comments, check the DL-DCL Matlab Simulation Video Master Code.%%%

clc
clear

T = 600; % discrete-time steps
deltaT = 0.1; % sampling rate

parfor k_N = 1:13
    N = k_N + 2; % total no. of robots in the system
  
    CLsf_pos = zeros(T,N,20);
    Tot_CLsf_pos = zeros(T,20);
    CLsf_pos_IMUf = zeros(T,N,20);
    Tot_CLsf_pos_IMUf = zeros(T,20);
    CLsf_pos_RSf = zeros(T,N,20);
    Tot_CLsf_pos_RSf = zeros(T,20);
    CLsf_pos_KF = zeros(T,N,20);
    Tot_CLsf_pos_KF = zeros(T,20);
    CLsf_pos_CI = zeros(T,N,20);
    Tot_CLsf_pos_CI = zeros(T,20);
    CLsf_pos_CU = zeros(T,N,20);
    Tot_CLsf_pos_CU = zeros(T,20);
    
for n = 1:20 % no. of simulation runs
%% Preliminary Parameters and variables

rng(n*10) % change seed for random no. generation

deltaT = 0.1; % sampling rate
xr = zeros(2,T,N); % position vectors of N robots over T time-steps
phir = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr = zeros(2,T,N);
delta_phir = zeros(1,T,N);

xr_gl = zeros(2,T,N);
phir_gl = zeros(1,T,N);

RotM = @(phi) [cos(phi) -sin(phi); sin(phi) cos(phi)]; % robot body-axis to global coordinate axis

create_containers = @(n)arrayfun(@(x)containers.Map(), 1:n, 'UniformOutput', false);

xr_IMUf = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_IMUf = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_IMUf = zeros(2,T,N);
delta_phir_IMUf = zeros(1,T,N);

xr_RSf = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_RSf = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_RSf = zeros(2,T,N);
delta_phir_RSf = zeros(1,T,N);

xr_KF = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_KF = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_KF = zeros(2,T,N);
delta_phir_KF = zeros(1,T,N);

xr_CI = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_CI = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_CI = zeros(2,T,N);
delta_phir_CI = zeros(1,T,N);

xr_CU = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_CU = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_CU = zeros(2,T,N);
delta_phir_CU = zeros(1,T,N);
%% Landmark Dynamics

xTg = zeros(2,T); % position vector of the Landmark over T time-steps
phiTg = zeros(1,T); % heading angle of the Landmark over T time-steps

wbar_Tg = zeros(1,T);
wbar_Tg(1,1:floor(T/6)) = 1.5*pi/(T*deltaT);
wbar_Tg(1,ceil(T/6):floor(2*T/6)) = 0;
wbar_Tg(1,ceil(2*T/6):floor(3*T/6)) = -1.5*pi/(T*deltaT);
wbar_Tg(1,ceil(3*T/6):floor(5*T/6)) = -3*pi/(T*deltaT);
wbar_Tg(1,ceil(5*T/6):ceil(6*T/6)) = -1*pi/(T*deltaT);

vbar_Tg = zeros(2,T);
vbar_Tg(:,1:floor(T/6)) = [linspace(0,1,length(1:floor(T/6))); linspace(0,0,length(1:floor(T/6)))];
vbar_Tg(:,ceil(T/6):floor(2*T/6)) = [linspace(0,0,length(ceil(T/6):floor(2*T/6))); linspace(-1,0,length(ceil(T/6):floor(2*T/6)))];
vbar_Tg(:,ceil(2*T/6):floor(3*T/6)) = [linspace(0,1,length(ceil(2*T/6):floor(3*T/6))); linspace(0,0,length(ceil(2*T/6):floor(3*T/6)))];
vbar_Tg(:,ceil(3*T/6):floor(5*T/6)) = [linspace(0,1,length(ceil(3*T/6):floor(5*T/6))); linspace(0,0,length(ceil(3*T/6):floor(5*T/6)))];
vbar_Tg(:,ceil(5*T/6):ceil(6*T/6)) = [linspace(1,0,length(ceil(5*T/6):ceil(6*T/6))); linspace(0,0,length(ceil(5*T/6):ceil(6*T/6)))];


xTg(:,1) = [10;10];
phiTg(:,1) = 0;
for t = 1:T-1
   [xTg(:,t+1),phiTg(:,t+1)] = targetctrldyn(deltaT,xTg(:,t),phiTg(:,t),RotM(phiTg(:,t)),vbar_Tg(:,t),wbar_Tg(:,t));
end

%% Comm. Network Dynamics

G = zeros(N);

for i=1:N
    for j=1:N
       if mod(abs(i-j),N-1) <= 1
          G(i,j) = 1; 
       else
          G(i,j) = 0; 
       end
    end 
end
 
% G = [1 1 1 0 1 1;
%      1 1 1 1 0 1;
%      1 1 1 1 1 0;
%      0 1 1 1 1 1;
%      1 0 1 1 1 1;
%      1 1 0 1 1 1];

%  
% figure()
% graph(G).plot

p = 0.0; % probability of range/comm. link failure
A = double(rand(N,N,T) > p);
A(:,:,1) = ones(N,N);
for k = 2:T
   A(:,:,k) = A(:,:,k) - diag(diag(A(:,:,k)));
   A(:,:,k) = A(:,:,k) - tril(A(:,:,k));
   A(:,:,k) = A(:,:,k) + A(:,:,k)' + diag(ones(1,N));
end

Gdyn = G.*A;

%% Robot dynamics
dTs_obrd = (2/6)*N;
for i = 1:N
   xr(:,1,i) = [10+dTs_obrd*cos(2*i*pi/N); 10+dTs_obrd*sin(2*i*pi/N)];
   xr_gl(:,1,i) = [10+dTs_obrd*cos(2*i*pi/N); 10+dTs_obrd*sin(2*i*pi/N)];
   xr_IMUf(:,1,i) = [10+dTs_obrd*cos(2*i*pi/N); 10+dTs_obrd*sin(2*i*pi/N)];
   xr_RSf(:,1,i) = [10+dTs_obrd*cos(2*i*pi/N); 10+dTs_obrd*sin(2*i*pi/N)];
   xr_KF(:,1,i) = [10+dTs_obrd*cos(2*i*pi/N); 10+dTs_obrd*sin(2*i*pi/N)];
   xr_CI(:,1,i) = [10+dTs_obrd*cos(2*i*pi/N); 10+dTs_obrd*sin(2*i*pi/N)];
   xr_CU(:,1,i) = [10+dTs_obrd*cos(2*i*pi/N); 10+dTs_obrd*sin(2*i*pi/N)];
   phir(:,1,i) = pi/2;
   phir_gl(:,1,i) = pi/2;
   phir_IMUf(:,1,i) = pi/2;
   phir_RSf(:,1,i) = pi/2;
   phir_KF(:,1,i) = pi/2;
   phir_CI(:,1,i) = pi/2;
   phir_CU(:,1,i) = pi/2;
end

dT_safe = (10*N/6)*ones(1,N);
kvTf = 6*ones(1,N);
kvNs = 12*ones(1,N);
kvRT = 0.5*ones(1,N);
kw = 10*ones(1,N);

nu_mean = zeros(T,N);
nu_mean(1:floor(T/6),:) = repmat([0.1 (1.1/2)*ones(1,floor(N/2)+mod(N,2)-1) 1 (1.1/2)*ones(1,floor(N/2)-1)],length(1:floor(T/6)),1);
nu_mean(ceil(T/6):ceil(6*T/6),:) = repmat([0.1 (5.1/2)*ones(1,floor(N/2)+mod(N,2)-1) 5 (5.1/2)*ones(1,floor(N/2)-1)],length(ceil(T/6):ceil(6*T/6)),1);
nu_var = [0.1 (4.1/2)*ones(1,floor(N/2)+mod(N,2)-1) 4 (4.1/2)*ones(1,floor(N/2)-1)];
xr_nu = nu_mean + repmat(sqrt(nu_var),T,1).*randn(T,N); % disturbance in the translational kinematics of robots

mu_mean = zeros(T,N);
mu_var = zeros(T,N);

mu_val01 = [0.1*ones(1,floor(N/2)+mod(N,2)-1) 3*ones(1,floor(N/2)+mod(N,2)-2)];
mu_val02 = [0.1*ones(1,floor(N/2)-1) 3*ones(1,floor(N/2)-1)];
mu_vec01ID = randperm(length(mu_val01),length(ones(1,floor(N/2)+mod(N,2)-1)))
mu_vec02ID = randperm(length(mu_val02),length(ones(1,floor(N/2)-1)))
mu_vec01 = mu_val01(mu_vec01ID)
mu_vec02 = mu_val02(mu_vec02ID)

mu_mean(1:floor(T/6),:) = repmat([0.1 (0.2/2)*ones(1,floor(N/2)+mod(N,2)-1) 0.1 (0.2/2)*ones(1,floor(N/2)-1)],length(1:floor(T/6)),1);
mu_mean(ceil(T/6):ceil(6*T/6),:) = repmat([0.1 mu_vec01.*ones(1,floor(N/2)+mod(N,2)-1) 6 mu_vec02.*ones(1,floor(N/2)-1)],length(ceil(T/6):ceil(6*T/6)),1);

mu_var(1:floor(T/6),:) = repmat([0.1 (1.1/2)*ones(1,floor(N/2)+mod(N,2)-1) 1 (1.1/2)*ones(1,floor(N/2)-1)],length(1:floor(T/6)),1);
mu_var(ceil(T/6):ceil(6*T/6),:) = repmat([0.1 mu_vec01.*ones(1,floor(N/2)+mod(N,2)-1) 6 mu_vec02.*ones(1,floor(N/2)-1)],length(ceil(T/6):ceil(6*T/6)),1);

mu_dx = 0.5*mu_mean + 0.05*sqrt(mu_var).*randn(T,N); % inter-robot relative distance noise in the range sensor

mu_cov = mu_mean.*(0.05)^2;
%
IMU_mean = zeros(T,N);
IMU_var = zeros(T,N);

IMU_val01 = [0.1*ones(1,floor(N/2)+mod(N,2)-1) 6*ones(1,floor(N/2)+mod(N,2)-2)];
IMU_val02 = [0.1*ones(1,floor(N/2)-1) 6*ones(1,floor(N/2)-1)];
IMU_vec01ID = randperm(length(IMU_val01),length(ones(1,floor(N/2)+mod(N,2)-1)))
IMU_vec02ID = randperm(length(IMU_val02),length(ones(1,floor(N/2)-1)))
IMU_vec01 = IMU_val01(IMU_vec01ID)
IMU_vec02 = IMU_val02(IMU_vec02ID)

IMU_mean(1:floor(T/6),:) = repmat([0.1 (0.2/2)*ones(1,floor(N/2)+mod(N,2)-1) 0.1 (0.2/2)*ones(1,floor(N/2)-1)],length(1:floor(T/6)),1);
IMU_mean(ceil(T/6):ceil(6*T/6),:) = repmat([0.1 3*ones(1,floor(N/2)+mod(N,2)-1) 6 3*ones(1,floor(N/2)-1)],length(ceil(T/6):ceil(6*T/6)),1);

IMU_var(1:floor(T/6),:) = repmat([0.1 (0.2/2)*ones(1,floor(N/2)+mod(N,2)-1) 0.1 (0.2/2)*ones(1,floor(N/2)-1)],length(1:floor(T/6)),1);
IMU_var(ceil(T/6):ceil(6*T/6),:) = repmat([0.1 3*ones(1,floor(N/2)+mod(N,2)-1) 6 3*ones(1,floor(N/2)-1)],length(ceil(T/6):ceil(6*T/6)),1);

IMU_noise = 0.5*IMU_mean + 0.05*sqrt(IMU_var).*randn(T,N); % noise in the IMU measurements

IMU_cov = IMU_mean.*(0.05)^2;
%% 
xr_net = zeros(2,T,N);
xr_IMU = zeros(2,T,N);
xr_Prj = zeros(2,T,N);
xr_Est = zeros(2,T,N);

phir_net = zeros(1,T,N);
phir_IMU = zeros(1,T,N);
phir_Prj = zeros(1,T,N);
phir_Est = zeros(1,T,N);

xr_net_IMUf = zeros(2,T,N);
phir_net_IMUf = zeros(1,T,N);

xr_net_RSf = zeros(2,T,N);
phir_net_RSf = zeros(1,T,N);

xr_net_KF = zeros(2,T,N);
phir_net_KF = zeros(1,T,N);
xr_IMU_KF = zeros(2,T,N);
phir_IMU_KF = zeros(1,T,N);

xr_net_CI = zeros(2,T,N);
phir_net_CI = zeros(1,T,N);
xr_IMU_CI = zeros(2,T,N);
phir_IMU_CI = zeros(1,T,N);

xr_net_CU = zeros(2,T,N);
phir_net_CU = zeros(1,T,N);
xr_IMU_CU = zeros(2,T,N);
phir_IMU_CU = zeros(1,T,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%L_Alpha = create_containers(N);
%L_Beta = create_containers(N);
L_Wij = create_containers(N);
L_Zi = zeros(T,N);

%Lpr_Alpha = create_containers(N);
%Lpr_Beta = create_containers(N);
Lpr_Wij = create_containers(N);
Lpr_Zi = zeros(T,N);


%L_Alpha_phi = create_containers(N);
%L_Beta_phi = create_containers(N);
L_Wij_phi = create_containers(N);
L_Zi_phi = zeros(T,N);

%Lpr_Alpha_phi = create_containers(N);
%Lpr_Beta_phi = create_containers(N);
Lpr_Wij_phi = create_containers(N);
Lpr_Zi_phi = zeros(T,N);

s_count = zeros(T,N);
s_count_phi = zeros(T,N);

for t=1:T
   
   for i = 1:N
          xr_IMU(:,t,i) = IMU_model(xr(:,t,i),IMU_mean(t,i));
          phir_IMU(:,t,i) = IMU_model_phi(phir(:,t,i),IMU_mean(t,i));
          if t > 1
              xr_Prj(:,t,i) = xr_Est(:,t-1,i) + delta_xr(:,t,i);
              phir_Prj(:,t,i) = wrapToPi(phir_Est(:,t-1,i) + delta_phir(:,t,i));
          else
              xr_Prj(:,t,i) = xr_IMU(:,t,i) + delta_xr(:,t,i);
              phir_Prj(:,t,i) = wrapToPi(phir_Est(:,t,i) + delta_phir(:,t,i));
          end
          
          xr_IMU_KF(:,t,i) = IMU_model(xr_KF(:,t,i),IMU_mean(t,i));
          phir_IMU_KF(:,t,i) = IMU_model_phi(phir_KF(:,t,i),IMU_mean(t,i));
          
          xr_IMU_CI(:,t,i) = IMU_model(xr_CI(:,t,i),IMU_mean(t,i));
          phir_IMU_CI(:,t,i) = IMU_model_phi(phir_CI(:,t,i),IMU_mean(t,i));
          
          xr_IMU_CU(:,t,i) = IMU_model(xr_CU(:,t,i),IMU_mean(t,i));
          phir_IMU_CU(:,t,i) = IMU_model_phi(phir_CU(:,t,i),IMU_mean(t,i));
   end 
    
   for i = 1:N
      %%%%%%%%%DL-DCL%%%%%%%%% 
      if t > 1 
      [xr_Est(:,t,i),s_count(t,i),L_Zi(t,i),Lpr_Zi(t,i)] = DLDCL1hop(t, s_count(t-1,i), N, Gdyn(:,:,t), i, permute(xr(:,t,:),[1 3 2]), permute(xr_Prj(:,t,:),[1 3 2]), permute(xr_IMU(:,t,:),[1 3 2]), xTg(:,t), mu_mean(t,:), L_Wij{i}, L_Zi(t-1,i), Lpr_Wij{i}, Lpr_Zi(t-1,i)); 
      else
      [xr_Est(:,t,i),s_count(t,i),L_Zi(t,i),Lpr_Zi(t,i)] = DLDCL1hop(t, 0, N, Gdyn(:,:,t), i, permute(xr(:,t,:),[1 3 2]), permute(xr_Prj(:,t,:),[1 3 2]), permute(xr_IMU(:,t,:),[1 3 2]), xTg(:,t), mu_mean(t,:), L_Wij{i}, 0, Lpr_Wij{i}, 0);     
      end
      
      %
      
      if t > 1 
      [phir_Est(:,t,i),s_count_phi(t,i),L_Zi_phi(t,i),Lpr_Zi_phi(t,i)] = DLDCL1hop_phi(t, s_count_phi(t-1,i), N, Gdyn(:,:,t), i, permute(phir(:,t,:),[1 3 2]), permute(phir_Prj(:,t,:),[1 3 2]), permute(phir_IMU(:,t,:),[1 3 2]), phiTg(:,t), mu_mean(t,:), L_Wij_phi{i}, L_Zi_phi(t-1,i), Lpr_Wij_phi{i}, Lpr_Zi_phi(t-1,i)); 
      else
      [phir_Est(:,t,i),s_count_phi(t,i),L_Zi_phi(t,i),Lpr_Zi_phi(t,i)] = DLDCL1hop_phi(t, 0, N, Gdyn(:,:,t), i, permute(phir(:,t,:),[1 3 2]), permute(phir_Prj(:,t,:),[1 3 2]), permute(phir_IMU(:,t,:),[1 3 2]), phiTg(:,t), mu_mean(t,:), L_Wij_phi{i}, 0, Lpr_Wij_phi{i}, 0);     
      end
      %%%%%%%%%DL-DCL%%%%%%%%%
      
      xr_net_KF(:,t,i) = DCLKF(N, Gdyn(:,:,t), i, permute(xr_KF(:,t,:),[1 3 2]), permute(xr_IMU_KF(:,t,:),[1 3 2]), xTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      phir_net_KF(:,t,i) = DCLKF_phi(N, Gdyn(:,:,t), i, permute(phir_KF(:,t,:),[1 3 2]), permute(phir_IMU_KF(:,t,:),[1 3 2]), phiTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      
      %%%
      
      xr_net_CI(:,t,i) = DCLCI(N, Gdyn(:,:,t), i, permute(xr_CI(:,t,:),[1 3 2]), permute(xr_IMU_CI(:,t,:),[1 3 2]), xTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      phir_net_CI(:,t,i) = DCLCI_phi(N, Gdyn(:,:,t), i, permute(phir_CI(:,t,:),[1 3 2]), permute(phir_IMU_CI(:,t,:),[1 3 2]), phiTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      
      %%%
      
      xr_net_CU(:,t,i) = DCLCU(N, Gdyn(:,:,t), i, permute(xr_CU(:,t,:),[1 3 2]), permute(xr_IMU_CU(:,t,:),[1 3 2]), xTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      phir_net_CU(:,t,i) = DCLCU_phi(N, Gdyn(:,:,t), i, permute(phir_CU(:,t,:),[1 3 2]), permute(phir_IMU_CU(:,t,:),[1 3 2]), phiTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      
   end
    
   %%if noise_set == 1 
        for i = 1:N
          xr_net_IMUf(:,t,i) = IMU_model(xr_IMUf(:,t,i),IMU_mean(t,i));  
          phir_net_IMUf(:,t,i) = IMU_model_phi(phir_IMUf(:,t,i),IMU_mean(t,i)); 
        end
   %%elseif noise_set == 2
        for i = 1:N
          xr_net_RSf(:,t,i) = xTg(:,t) - RangeSensing(xr_RSf(:,t,i),xTg(:,t),mu_mean(t,i));  
          phir_net_RSf(:,t,i) = phiTg(:,t) - RelPhiSens(phir_RSf(:,t,i),phiTg(:,t),mu_mean(t,i));
        end  
   %%elseif noise_set == 3
          xr_net(:,t,:) = xr_Est(:,t,:); 
          phir_net(:,t,:) = phir_Est(:,t,:); 
   %%else
%           xr_net(:,t,:) = xr(:,t,:);  
%           phir_net(:,t,:) = phir(:,t,:); 
   %%end
      
   i = 1;
   ID_prev = N;
   ID_nxt = i+1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr(:,t,:),[1 3 2]));
   [xr(:,t+1,i),phir(:,t+1,i),delta_xr(:,t+1,i),delta_phir(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i),kvRT(i), kw(i), deltaT, xr_net(:,t,i), phir_net(:,t,i), RotM(phir(:,t,i)), xr_net(:,t,ID_prev), xr_net(:,t,ID_nxt),xr_net(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr(:,t,i),phir(:,t,i));
   for i = 2:N-1
       ID_prev = i-1;
       ID_nxt = i+1;
       [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr(:,t,:),[1 3 2]));
       [xr(:,t+1,i),phir(:,t+1,i),delta_xr(:,t+1,i),delta_phir(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net(:,t,i), phir_net(:,t,i), RotM(phir(:,t,i)), xr_net(:,t,ID_prev), xr_net(:,t,ID_nxt),xr_net(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr(:,t,i),phir(:,t,i));
   end
   i = N;
   ID_prev = i-1;
   ID_nxt = 1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr(:,t,:),[1 3 2]));
   [xr(:,t+1,i),phir(:,t+1,i),delta_xr(:,t+1,i),delta_phir(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net(:,t,i), phir_net(:,t,i), RotM(phir(:,t,i)), xr_net(:,t,ID_prev), xr_net(:,t,ID_nxt),xr_net(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr(:,t,i),phir(:,t,i));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%
   
   i = 1;
   ID_prev = N;
   ID_nxt = i+1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_IMUf(:,t,:),[1 3 2]));
   [xr_IMUf(:,t+1,i),phir_IMUf(:,t+1,i),delta_xr_IMUf(:,t+1,i),delta_phir_IMUf(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i),kvRT(i), kw(i), deltaT, xr_net_IMUf(:,t,i), phir_net_IMUf(:,t,i), RotM(phir_IMUf(:,t,i)), xr_net_IMUf(:,t,ID_prev), xr_net_IMUf(:,t,ID_nxt),xr_net_IMUf(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_IMUf(:,t,i),phir_IMUf(:,t,i));
   for i = 2:N-1
       ID_prev = i-1;
       ID_nxt = i+1;
       [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_IMUf(:,t,:),[1 3 2]));
       [xr_IMUf(:,t+1,i),phir_IMUf(:,t+1,i),delta_xr_IMUf(:,t+1,i),delta_phir_IMUf(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_IMUf(:,t,i), phir_net_IMUf(:,t,i), RotM(phir_IMUf(:,t,i)), xr_net_IMUf(:,t,ID_prev), xr_net_IMUf(:,t,ID_nxt),xr_net_IMUf(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_IMUf(:,t,i),phir_IMUf(:,t,i));
   end
   i = N;
   ID_prev = i-1;
   ID_nxt = 1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_IMUf(:,t,:),[1 3 2]));
   [xr_IMUf(:,t+1,i),phir_IMUf(:,t+1,i),delta_xr_IMUf(:,t+1,i),delta_phir_IMUf(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_IMUf(:,t,i), phir_net_IMUf(:,t,i), RotM(phir_IMUf(:,t,i)), xr_net_IMUf(:,t,ID_prev), xr_net_IMUf(:,t,ID_nxt),xr_net_IMUf(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_IMUf(:,t,i),phir_IMUf(:,t,i));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   i = 1;
   ID_prev = N;
   ID_nxt = i+1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_RSf(:,t,:),[1 3 2]));
   [xr_RSf(:,t+1,i),phir_RSf(:,t+1,i),delta_xr_RSf(:,t+1,i),delta_phir_RSf(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i),kvRT(i), kw(i), deltaT, xr_net_RSf(:,t,i), phir_net_RSf(:,t,i), RotM(phir_RSf(:,t,i)), xr_net_RSf(:,t,ID_prev), xr_net_RSf(:,t,ID_nxt),xr_net_RSf(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_RSf(:,t,i),phir_RSf(:,t,i));
   for i = 2:N-1
       ID_prev = i-1;
       ID_nxt = i+1;
       [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_RSf(:,t,:),[1 3 2]));
       [xr_RSf(:,t+1,i),phir_RSf(:,t+1,i),delta_xr_RSf(:,t+1,i),delta_phir_RSf(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_RSf(:,t,i), phir_net_RSf(:,t,i), RotM(phir_RSf(:,t,i)), xr_net_RSf(:,t,ID_prev), xr_net_RSf(:,t,ID_nxt),xr_net_RSf(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_RSf(:,t,i),phir_RSf(:,t,i));
   end
   i = N;
   ID_prev = i-1;
   ID_nxt = 1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_RSf(:,t,:),[1 3 2]));
   [xr_RSf(:,t+1,i),phir_RSf(:,t+1,i),delta_xr_RSf(:,t+1,i),delta_phir_RSf(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_RSf(:,t,i), phir_net_RSf(:,t,i), RotM(phir_RSf(:,t,i)), xr_net_RSf(:,t,ID_prev), xr_net_RSf(:,t,ID_nxt),xr_net_RSf(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_RSf(:,t,i),phir_RSf(:,t,i));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   i = 1;
   ID_prev = N;
   ID_nxt = i+1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_KF(:,t,:),[1 3 2]));
   [xr_KF(:,t+1,i),phir_KF(:,t+1,i),delta_xr_KF(:,t+1,i),delta_phir_KF(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i),kvRT(i), kw(i), deltaT, xr_net_KF(:,t,i), phir_net_KF(:,t,i), RotM(phir_KF(:,t,i)), xr_net_KF(:,t,ID_prev), xr_net_KF(:,t,ID_nxt),xr_net_KF(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_KF(:,t,i),phir_KF(:,t,i));
   for i = 2:N-1
       ID_prev = i-1;
       ID_nxt = i+1;
       [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_KF(:,t,:),[1 3 2]));
       [xr_KF(:,t+1,i),phir_KF(:,t+1,i),delta_xr_KF(:,t+1,i),delta_phir_KF(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_KF(:,t,i), phir_net_KF(:,t,i), RotM(phir_KF(:,t,i)), xr_net_KF(:,t,ID_prev), xr_net_KF(:,t,ID_nxt),xr_net_KF(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_KF(:,t,i),phir_KF(:,t,i));
   end
   i = N;
   ID_prev = i-1;
   ID_nxt = 1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_KF(:,t,:),[1 3 2]));
   [xr_KF(:,t+1,i),phir_KF(:,t+1,i),delta_xr_KF(:,t+1,i),delta_phir_KF(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_KF(:,t,i), phir_net_KF(:,t,i), RotM(phir_KF(:,t,i)), xr_net_KF(:,t,ID_prev), xr_net_KF(:,t,ID_nxt),xr_net_KF(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_KF(:,t,i),phir_KF(:,t,i));
end


%%%

i = 1;
   ID_prev = N;
   ID_nxt = i+1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_CI(:,t,:),[1 3 2]));
   [xr_CI(:,t+1,i),phir_CI(:,t+1,i),delta_xr_CI(:,t+1,i),delta_phir_CI(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i),kvRT(i), kw(i), deltaT, xr_net_CI(:,t,i), phir_net_CI(:,t,i), RotM(phir_CI(:,t,i)), xr_net_CI(:,t,ID_prev), xr_net_CI(:,t,ID_nxt),xr_net_CI(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_CI(:,t,i),phir_CI(:,t,i));
   for i = 2:N-1
       ID_prev = i-1;
       ID_nxt = i+1;
       [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_CI(:,t,:),[1 3 2]));
       [xr_CI(:,t+1,i),phir_CI(:,t+1,i),delta_xr_CI(:,t+1,i),delta_phir_CI(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_CI(:,t,i), phir_net_CI(:,t,i), RotM(phir_CI(:,t,i)), xr_net_CI(:,t,ID_prev), xr_net_CI(:,t,ID_nxt),xr_net_CI(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_CI(:,t,i),phir_CI(:,t,i));
   end
   i = N;
   ID_prev = i-1;
   ID_nxt = 1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_CI(:,t,:),[1 3 2]));
   [xr_CI(:,t+1,i),phir_CI(:,t+1,i),delta_xr_CI(:,t+1,i),delta_phir_CI(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_CI(:,t,i), phir_net_CI(:,t,i), RotM(phir_CI(:,t,i)), xr_net_CI(:,t,ID_prev), xr_net_CI(:,t,ID_nxt),xr_net_CI(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_CI(:,t,i),phir_CI(:,t,i));
   
   %%%
   
   i = 1;
   ID_prev = N;
   ID_nxt = i+1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_CU(:,t,:),[1 3 2]));
   [xr_CU(:,t+1,i),phir_CU(:,t+1,i),delta_xr_CU(:,t+1,i),delta_phir_CU(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i),kvRT(i), kw(i), deltaT, xr_net_CU(:,t,i), phir_net_CU(:,t,i), RotM(phir_CU(:,t,i)), xr_net_CU(:,t,ID_prev), xr_net_CU(:,t,ID_nxt),xr_net_CU(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_CU(:,t,i),phir_CU(:,t,i));
   for i = 2:N-1
       ID_prev = i-1;
       ID_nxt = i+1;
       [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_CU(:,t,:),[1 3 2]));
       [xr_CU(:,t+1,i),phir_CU(:,t+1,i),delta_xr_CU(:,t+1,i),delta_phir_CU(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_CU(:,t,i), phir_net_CU(:,t,i), RotM(phir_CU(:,t,i)), xr_net_CU(:,t,ID_prev), xr_net_CU(:,t,ID_nxt),xr_net_CU(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_CU(:,t,i),phir_CU(:,t,i));
   end
   i = N;
   ID_prev = i-1;
   ID_nxt = 1;
   [ID_nrst_rtb,~] = NrstNeighDist(i,permute(xr_CU(:,t,:),[1 3 2]));
   [xr_CU(:,t+1,i),phir_CU(:,t+1,i),delta_xr_CU(:,t+1,i),delta_phir_CU(:,t+1,i)] = robotctrldyNnoise(dT_safe(i), kvTf(i), kvNs(i), kvRT(i), kw(i), deltaT, xr_net_CU(:,t,i), phir_net_CU(:,t,i), RotM(phir_CU(:,t,i)), xr_net_CU(:,t,ID_prev), xr_net_CU(:,t,ID_nxt),xr_net_CU(:,t,ID_nrst_rtb), xTg(:,t), xr_nu(t,i),ID_prev,ID_nxt,ID_nrst_rtb,xr_CU(:,t,i),phir_CU(:,t,i));

%%

lssfn = @(x,y) min(sqrt(sum((x-y).^2))./15,1);

lss_pos = zeros(T,N);
Clss_pos = zeros(T,N);

for t = 1:T
   lss_pos(t,:) = lssfn(permute(xr(:,t,:),[1 3 2]),permute(xr_net(:,t,:),[1 3 2])); 
   if t > 1
        Clss_pos(t,:) = Clss_pos(t-1,:) + lss_pos(t,:);  
   else
        Clss_pos(t,:) = 0 + lss_pos(t,:); 
   end
end

Tot_Clss_pos = sum(Clss_pos,2);

CLsf_pos(:,:,n) = Clss_pos;
Tot_CLsf_pos(:,n) = Tot_Clss_pos;

%%%%%%%%%%%%%%%%%%%%%%%%

lss_pos_IMUf = zeros(T,N);
Clss_pos_IMUf = zeros(T,N);

for t = 1:T
   lss_pos_IMUf(t,:) = lssfn(permute(xr_IMUf(:,t,:),[1 3 2]),permute(xr_net_IMUf(:,t,:),[1 3 2])); 
   if t > 1
        Clss_pos_IMUf(t,:) = Clss_pos_IMUf(t-1,:) + lss_pos_IMUf(t,:);  
   else
        Clss_pos_IMUf(t,:) = 0 + lss_pos_IMUf(t,:); 
   end
end

Tot_Clss_pos_IMUf = sum(Clss_pos_IMUf,2);

CLsf_pos_IMUf(:,:,n) = Clss_pos_IMUf;
Tot_CLsf_pos_IMUf(:,n) = Tot_Clss_pos_IMUf;

%%%%%%%%%%%%%%%%%%%%%%%%%

lss_pos_RSf = zeros(T,N);
Clss_pos_RSf = zeros(T,N);

for t = 1:T
   lss_pos_RSf(t,:) = lssfn(permute(xr_RSf(:,t,:),[1 3 2]),permute(xr_net_RSf(:,t,:),[1 3 2])); 
   if t > 1
        Clss_pos_RSf(t,:) = Clss_pos_RSf(t-1,:) + lss_pos_RSf(t,:);  
   else
        Clss_pos_RSf(t,:) = 0 + lss_pos_RSf(t,:); 
   end
end

Tot_Clss_pos_RSf = sum(Clss_pos_RSf,2);

CLsf_pos_RSf(:,:,n) = Clss_pos_RSf;
Tot_CLsf_pos_RSf(:,n) = Tot_Clss_pos_RSf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lss_pos_KF = zeros(T,N);
Clss_pos_KF = zeros(T,N);

for t = 1:T
   lss_pos_KF(t,:) = lssfn(permute(xr_KF(:,t,:),[1 3 2]),permute(xr_net_KF(:,t,:),[1 3 2])); 
   if t > 1
        Clss_pos_KF(t,:) = Clss_pos_KF(t-1,:) + lss_pos_KF(t,:);  
   else
        Clss_pos_KF(t,:) = 0 + lss_pos_KF(t,:); 
   end
end

Tot_Clss_pos_KF = sum(Clss_pos_KF,2);

CLsf_pos_KF(:,:,n) = Clss_pos_KF;
Tot_CLsf_pos_KF(:,n) = Tot_Clss_pos_KF;

%%%%%%%%%%%%%%%%%%%%%%%

lss_pos_CI = zeros(T,N);
Clss_pos_CI = zeros(T,N);

for t = 1:T
   lss_pos_CI(t,:) = lssfn(permute(xr_CI(:,t,:),[1 3 2]),permute(xr_net_CI(:,t,:),[1 3 2])); 
   if t > 1
        Clss_pos_CI(t,:) = Clss_pos_CI(t-1,:) + lss_pos_CI(t,:);  
   else
        Clss_pos_CI(t,:) = 0 + lss_pos_CI(t,:); 
   end
end

Tot_Clss_pos_CI = sum(Clss_pos_CI,2);

CLsf_pos_CI(:,:,n) = Clss_pos_CI;
Tot_CLsf_pos_CI(:,n) = Tot_Clss_pos_CI;


%%%%%%%%%%%%%%%%%%%%%%%

lss_pos_CU = zeros(T,N);
Clss_pos_CU = zeros(T,N);

for t = 1:T
   lss_pos_CU(t,:) = lssfn(permute(xr_CU(:,t,:),[1 3 2]),permute(xr_net_CU(:,t,:),[1 3 2])); 
   if t > 1
        Clss_pos_CU(t,:) = Clss_pos_CU(t-1,:) + lss_pos_CU(t,:);  
   else
        Clss_pos_CU(t,:) = 0 + lss_pos_CU(t,:); 
   end
end

Tot_Clss_pos_CU = sum(Clss_pos_CU,2);

CLsf_pos_CU(:,:,n) = Clss_pos_CU;
Tot_CLsf_pos_CU(:,n) = Tot_Clss_pos_CU;

end
%%

nrbt = ceil(N/2) + 1;

CL_pos = permute(CLsf_pos,[1 3 2]);
CL_pos_nrbt = mean(CL_pos(:,:,nrbt),2);

Tot_CL_pos = mean(Tot_CLsf_pos,2);

%%%%%%%%%%%%%%

CL_pos_IMUf = permute(CLsf_pos_IMUf,[1 3 2]);
CL_pos_nrbt_IMUf = mean(CL_pos_IMUf(:,:,nrbt),2);

Tot_CL_pos_IMUf = mean(Tot_CLsf_pos_IMUf,2);

%%%%%%%%%%%%%%

CL_pos_RSf = permute(CLsf_pos_RSf,[1 3 2]);
CL_pos_nrbt_RSf = mean(CL_pos_RSf(:,:,nrbt),2);

Tot_CL_pos_RSf = mean(Tot_CLsf_pos_RSf,2);

%%%%%%%%%%%%%%

CL_pos_KF = permute(CLsf_pos_KF,[1 3 2]);
CL_pos_nrbt_KF = mean(CL_pos_KF(:,:,nrbt),2);

Tot_CL_pos_KF = mean(Tot_CLsf_pos_KF,2);

%%%%%%%%%%%%%%

CL_pos_CI = permute(CLsf_pos_CI,[1 3 2]);
CL_pos_nrbt_CI = mean(CL_pos_CI(:,:,nrbt),2);

Tot_CL_pos_CI = mean(Tot_CLsf_pos_CI,2);

%%%%%%%%%%%%%%

CL_pos_CU = permute(CLsf_pos_CU,[1 3 2]);
CL_pos_nrbt_CU = mean(CL_pos_CU(:,:,nrbt),2);

Tot_CL_pos_CU = mean(Tot_CLsf_pos_CU,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Avg_CL_pos(k_N) = Tot_CL_pos(end)/N;
Avg_CL_pos_IMUf(k_N) = Tot_CL_pos_IMUf(end)/N;
Avg_CL_pos_RSf(k_N) = Tot_CL_pos_RSf(end)/N;
Avg_CL_pos_KF(k_N) = Tot_CL_pos_KF(end)/N;
Avg_CL_pos_CI(k_N) = Tot_CL_pos_CI(end)/N;
Avg_CL_pos_CU(k_N) = Tot_CL_pos_CU(end)/N;

end
%%

figure()
plot(3:15,Avg_CL_pos,3:15,Avg_CL_pos_IMUf,3:15,Avg_CL_pos_RSf,3:15,Avg_CL_pos_KF,3:15,Avg_CL_pos_CI,3:15,Avg_CL_pos_CU,'LineWidth',1.2)
set(gca,'FontSize',12)
xlabel('Total no. of robots (N)','FontSize',12,'Interpreter','Latex');
ylabel('Avg. Cumul. Estim. Loss of all robots after 60 sec.','FontSize',12,'Interpreter','Latex');
title("Avg. cumul. estim. loss of all robots after 60 sec. vs total no. of robots",'FontSize',12,'Interpreter','Latex');
legend('DL-DCL','No Comm., IMU only','No Comm., RPSS only','KF','CI','CU','Interpreter','Latex')
grid on
grid minor

figure()
subplot(2,1,1)
plot(3:15,Avg_CL_pos,3:15,Avg_CL_pos_KF,3:15,Avg_CL_pos_CI,3:15,Avg_CL_pos_CU,'LineWidth',1.2)
set(gca,'FontSize',12)
xlabel('Total no. of robots (N)','FontSize',12,'Interpreter','Latex');
ylabel('Avg. Cumul. Estim. Loss','FontSize',12,'Interpreter','Latex');
title("Avg. cumul. estim. loss of all robots after 60 sec. vs total no. of robots",'FontSize',12,'Interpreter','Latex');
legend('DL-DCL','KF','CI','CU','Interpreter','Latex')
grid on
xlim([3 15])
grid minor

subplot(2,1,2)
plot([1:20], 1./[1:20],'LineWidth',1.2)
set(gca,'FontSize',12)
xlabel('Total No. of Robots: $N$','FontSize',12,'Interpreter','Latex');
ylabel('Reliability Cost','FontSize',12,'Interpreter','Latex');
title("Reliability Cost vs Total no. of Robots",'FontSize',12,'Interpreter','Latex');
grid on
grid minor
xlim([3 15])
grid minor