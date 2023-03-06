%%% DL-DCL Matlab Simulation Video Master Code; comparative study with KF, CI, CU %%%
%%% by Shubhankar Gupta, Suresh Sundaram %%%
%%% AIRL, Dept. of Aerospace Engg., IISc. Bangalore, India %%%
%%% For more details, please check the associated AAAI-23 conference paper titled,
%%% 'Moving-Landmark assisted Distributed Learning based Decentralized
%%% Cooperative Localization (DL-DCL) with Fault Tolerance.' %%%

clc
clear

n = 1;
%% Preliminary parameters and variables 

T = 1400; % discrete-time steps
N = 6; % total no. of robots
deltaT = 0.1; % sampling rate

%pose estimates from DL-DCL
xr = zeros(2,T,N); % position vectors of N robots over T time-steps
phir = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr = zeros(2,T,N);
delta_phir = zeros(1,T,N);

RotM = @(phi) [cos(phi) -sin(phi); sin(phi) cos(phi)]; % robot body-axis to global coordinate axis

create_containers = @(n)arrayfun(@(x)containers.Map(), 1:n, 'UniformOutput', false);

% pose esimates from IMU only
xr_IMUf = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_IMUf = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_IMUf = zeros(2,T,N);
delta_phir_IMUf = zeros(1,T,N);

%pose estimates from Relative Pose Sensing only
xr_RSf = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_RSf = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_RSf = zeros(2,T,N);
delta_phir_RSf = zeros(1,T,N);

%pose estimates from Kalman Fusion
xr_KF = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_KF = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_KF = zeros(2,T,N);
delta_phir_KF = zeros(1,T,N);

%pose estimates from Covariance Intersection
xr_CI = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_CI = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_CI = zeros(2,T,N);
delta_phir_CI = zeros(1,T,N);

%pose estimates from Covariance Union
xr_CU = zeros(2,T,N); % position vectors of N robots over T time-steps
phir_CU = zeros(1,T,N); % heading angle of N robots over T time-steps
delta_xr_CU = zeros(2,T,N);
delta_phir_CU = zeros(1,T,N);

%% Landmark Dynamics

xTg = zeros(2,T); % position vector of the Landmark over T time-steps
phiTg = zeros(1,T); % heading angle of the Landmark over T time-steps

%yaw-rate control sequence
wbar_Tg = zeros(1,T);
wbar_Tg(1,1:floor(T/6)) = 3*pi/(T*deltaT);
wbar_Tg(1,ceil(T/6):floor(2*T/6)) = 0;
wbar_Tg(1,ceil(2*T/6):floor(3*T/6)) = -3*pi/(T*deltaT);
wbar_Tg(1,ceil(3*T/6):floor(5*T/6)) = -6*pi/(T*deltaT);
wbar_Tg(1,ceil(5*T/6):ceil(6*T/6)) = -2*pi/(T*deltaT);

%velocity control sequence
vbar_Tg = zeros(2,T);
vbar_Tg(:,1:floor(T/6)) = [linspace(0,2,length(1:floor(T/6))); linspace(0,0,length(1:floor(T/6)))];
vbar_Tg(:,ceil(T/6):floor(2*T/6)) = [linspace(0,0,length(ceil(T/6):floor(2*T/6))); linspace(-2,0,length(ceil(T/6):floor(2*T/6)))];
vbar_Tg(:,ceil(2*T/6):floor(3*T/6)) = [linspace(0,2,length(ceil(2*T/6):floor(3*T/6))); linspace(0,0,length(ceil(2*T/6):floor(3*T/6)))];
vbar_Tg(:,ceil(3*T/6):floor(5*T/6)) = [linspace(0,2,length(ceil(3*T/6):floor(5*T/6))); linspace(0,0,length(ceil(3*T/6):floor(5*T/6)))];
vbar_Tg(:,ceil(5*T/6):ceil(6*T/6)) = [linspace(2,0,length(ceil(5*T/6):ceil(6*T/6))); linspace(0,0,length(ceil(5*T/6):ceil(6*T/6)))];

%landmark pose sequence 
xTg(:,1) = [10;10];
phiTg(:,1) = 0;
for t = 1:T-1
   [xTg(:,t+1),phiTg(:,t+1)] = targetctrldyn(deltaT,xTg(:,t),phiTg(:,t),RotM(phiTg(:,t)),vbar_Tg(:,t),wbar_Tg(:,t));
end

%% Interaction ('communication + observation' link) Network Dynamics

% base network topology for the random link-drop interaction network
G = [1 1 0 0 0 1;
     1 1 1 0 0 0;
     0 1 1 1 0 0;
     0 0 1 1 1 0;
     0 0 0 1 1 1;
     1 0 0 0 1 1];
 
p = 0.5; % probability of an interaction link failure
A = double(rand(N,N,T) > p);
A(:,:,1) = ones(N,N);
for k = 2:T
   A(:,:,k) = A(:,:,k) - diag(diag(A(:,:,k)));
   A(:,:,k) = A(:,:,k) - tril(A(:,:,k));
   A(:,:,k) = A(:,:,k) + A(:,:,k)' + diag(ones(1,N));
end

Gdyn = G.*A; % dynamic network sequence for the random interaction network

%% Robot dynamics

% initialization for the pose of the robots 
for i = 1:N
   xr(:,1,i) = [10+2*cos(2*i*pi/N); 10+2*sin(2*i*pi/N)];
   xr_IMUf(:,1,i) = [10+2*cos(2*i*pi/N); 10+2*sin(2*i*pi/N)];
   xr_RSf(:,1,i) = [10+2*cos(2*i*pi/N); 10+2*sin(2*i*pi/N)];
   xr_KF(:,1,i) = [10+2*cos(2*i*pi/N); 10+2*sin(2*i*pi/N)];
   xr_CI(:,1,i) = [10+2*cos(2*i*pi/N); 10+2*sin(2*i*pi/N)];
   xr_CU(:,1,i) = [10+2*cos(2*i*pi/N); 10+2*sin(2*i*pi/N)];
   phir(:,1,i) = pi/2;
   phir_IMUf(:,1,i) = pi/2;
   phir_RSf(:,1,i) = pi/2;
   phir_KF(:,1,i) = pi/2;
   phir_CI(:,1,i) = pi/2;
   phir_CU(:,1,i) = pi/2;
end

% control parameters
dT_safe = 10*ones(1,N);
kvTf = 6*ones(1,N);
kvNs = 12*ones(1,N);
kvRT = 0.5*ones(1,N);
kw = 10*ones(1,N);

% disturbances (mean and covariance) in the translational kinematics of robots
nu_mean = zeros(T,N);
nu_mean(1:floor(T/6),:) = repmat([0.1 1 1 2 4 5],length(1:floor(T/6)),1);
nu_mean(ceil(T/6):floor(2*T/6),:) = repmat([1 1 2 1 3 6],length(ceil(T/6):floor(2*T/6)),1);
nu_mean(ceil(2*T/6):floor(3*T/6),:) = repmat([3 3 4 5 3 3],length(ceil(2*T/6):floor(3*T/6)),1);
nu_mean(ceil(3*T/6):floor(5*T/6),:) = repmat([6 3 2 2 3 0.1],length(ceil(3*T/6):floor(5*T/6)),1);
nu_mean(ceil(5*T/6):ceil(6*T/6),:) = repmat([5 3 2 2 1 0.1],length(ceil(5*T/6):ceil(6*T/6)),1);
nu_var = [4 6 9 2 1 4];
xr_nu = nu_mean + repmat(sqrt(nu_var),T,1).*randn(T,N); 

% inter-robot relative pose sensing system (RPSS) noise (mean and covariance) 
mu_mean = zeros(T,N);
mu_var = zeros(T,N);

mu_val = [0.1 0.1 0.1 0.1 0.1 4 4 4];

mu_mean(1:floor(T/6),:) = repmat(mu_val(randperm(length(mu_val),6)),length(1:floor(T/6)),1);
mu_mean(ceil(T/6):floor(2*T/6),:) = repmat(mu_val(randperm(length(mu_val),6)),length(ceil(T/6):floor(2*T/6)),1);
mu_mean(ceil(2*T/6):floor(3*T/6),:) = repmat(mu_val(randperm(length(mu_val),6)),length(ceil(2*T/6):floor(3*T/6)),1);
mu_mean(ceil(3*T/6):floor(5*T/6),:) = repmat(mu_val(randperm(length(mu_val),6)),length(ceil(3*T/6):floor(5*T/6)),1);
mu_mean(ceil(5*T/6):ceil(6*T/6),:) = repmat(mu_val(randperm(length(mu_val),6)),length(ceil(5*T/6):ceil(6*T/6)),1);

mu_var(1:floor(T/6),:) = repmat([0.1 0.1 4 0.1 0.1 4],length(1:floor(T/6)),1);
mu_var(ceil(T/6):floor(2*T/6),:) = repmat([0.1 0.1 4 4 0.1 4],length(ceil(T/6):floor(2*T/6)),1);
mu_var(ceil(2*T/6):floor(3*T/6),:) = repmat([0.1 0.1 0.1 4 0.1 4],length(ceil(2*T/6):floor(3*T/6)),1);
mu_var(ceil(3*T/6):floor(5*T/6),:) = repmat([0.1 4 0.1 4 0.1 4],length(ceil(3*T/6):floor(5*T/6)),1);
mu_var(ceil(5*T/6):ceil(6*T/6),:) = repmat([0.1 4 0.1 4 0.1 4],length(ceil(5*T/6):ceil(6*T/6)),1);

mu_dx = 0.5*mu_mean + 0.05*sqrt(mu_var).*randn(T,N); 

mu_cov = mu_mean.*(0.05)^2;

% noise (mean and covariance) in the IMU measurements
IMU_mean = zeros(T,N);
IMU_var = zeros(T,N);

IMU_ord = randperm(6,3);

IMU_mean(1:floor(T/6),:) = repmat([0.1 0.1 0.1 0.1 0.1 0.1],length(1:floor(T/6)),1);
IMU_mean(ceil(T/6):floor(2*T/6),:) = repmat([0.1 0.1 0.1 0.1 0.1 0.1],length(ceil(T/6):floor(2*T/6)),1);
IMU_mean(ceil(2*T/6):floor(3*T/6),:) = repmat([0.1 0.1 0.1 0.1 0.1 0.1],length(ceil(2*T/6):floor(3*T/6)),1);
IMU_mean(ceil(3*T/6):floor(5*T/6),:) = repmat([0.1 0.1 0.1 0.1 0.1 0.1],length(ceil(3*T/6):floor(5*T/6)),1);
IMU_mean(ceil(5*T/6):ceil(6*T/6),:) = repmat([0.1 0.1 0.1 0.1 0.1 0.1],length(ceil(5*T/6):ceil(6*T/6)),1);

IMU_mean(1:T,IMU_ord(1)) = repmat(6,length(1:T),1);
IMU_mean(ceil(T/6):T,IMU_ord(2)) = repmat(6,length(ceil(T/6):T),1);
IMU_mean(ceil(2*T/6):T,IMU_ord(3)) = repmat(6,length(ceil(2*T/6):T),1);

IMU_var(1:floor(T/6),:) = repmat([0.1 0.1 0.1 0.1 0.1 4],length(1:floor(T/6)),1);
IMU_var(ceil(T/6):floor(2*T/6),:) = repmat([0.1 0.1 0.1 0.1 4 5],length(ceil(T/6):floor(2*T/6)),1);
IMU_var(ceil(2*T/6):floor(3*T/6),:) = repmat([0.1 0.1 0.1 0.5 5 6],length(ceil(2*T/6):floor(3*T/6)),1);
IMU_var(ceil(3*T/6):floor(5*T/6),:) = repmat([0.1 0.1 0.1 4 5 6],length(ceil(3*T/6):floor(5*T/6)),1);
IMU_var(ceil(5*T/6):ceil(6*T/6),:) = repmat([0.1 0.1 0.1 5 6 6],length(ceil(5*T/6):ceil(6*T/6)),1);

IMU_noise = 0.5*IMU_mean + 0.05*sqrt(IMU_var).*randn(T,N); 

IMU_cov = IMU_mean.*(0.05)^2;

%% 

% pose sequence variables decalaration
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

%%% containers for storing cumulative estimation loss values %%%%%%%%%%%%%%
L_Wij = create_containers(N);
L_Zi = zeros(T,N);


Lpr_Wij = create_containers(N);
Lpr_Zi = zeros(T,N);


L_Wij_phi = create_containers(N);
L_Zi_phi = zeros(T,N);


Lpr_Wij_phi = create_containers(N);
Lpr_Zi_phi = zeros(T,N);

s_count = zeros(T,N);
s_count_phi = zeros(T,N);

%%% Simulation Loop
for t=1:T
   
   for i = 1:N
          % get the IMU estimates (x_IMU) (DL-DCL)
          xr_IMU(:,t,i) = IMU_model(xr(:,t,i),IMU_mean(t,i)); 
          phir_IMU(:,t,i) = IMU_model_phi(phir(:,t,i),IMU_mean(t,i));
          
          % get the projected estimates (x_PRJ) (DL-DCL)
          if t > 1
              xr_Prj(:,t,i) = xr_Est(:,t-1,i) + delta_xr(:,t,i);
              phir_Prj(:,t,i) = wrapToPi(phir_Est(:,t-1,i) + delta_phir(:,t,i));
          else
              xr_Prj(:,t,i) = xr_IMU(:,t,i) + delta_xr(:,t,i);
              phir_Prj(:,t,i) = wrapToPi(phir_Est(:,t,i) + delta_phir(:,t,i));
          end
          
          % estimates for simulation of other methods: KF, CI, CU
          xr_IMU_KF(:,t,i) = IMU_model(xr_KF(:,t,i),IMU_mean(t,i));
          phir_IMU_KF(:,t,i) = IMU_model_phi(phir_KF(:,t,i),IMU_mean(t,i));
          
          xr_IMU_CI(:,t,i) = IMU_model(xr_CI(:,t,i),IMU_mean(t,i));
          phir_IMU_CI(:,t,i) = IMU_model_phi(phir_CI(:,t,i),IMU_mean(t,i));
          
          xr_IMU_CU(:,t,i) = IMU_model(xr_CU(:,t,i),IMU_mean(t,i));
          phir_IMU_CU(:,t,i) = IMU_model_phi(phir_CU(:,t,i),IMU_mean(t,i));
   end 
   
   % run each fusion method for each robot
   for i = 1:N
      %%%%%%%%%DL-DCL%%%%%%%%%%%%%%%%%%%%% 
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
      
      %%%%%%%%%Kalman Fusion%%%%%%%%%%%%%%
      xr_net_KF(:,t,i) = DCLKF(N, Gdyn(:,:,t), i, permute(xr_KF(:,t,:),[1 3 2]), permute(xr_IMU_KF(:,t,:),[1 3 2]), xTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      phir_net_KF(:,t,i) = DCLKF_phi(N, Gdyn(:,:,t), i, permute(phir_KF(:,t,:),[1 3 2]), permute(phir_IMU_KF(:,t,:),[1 3 2]), phiTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      
      
      %%%%%%%%%Covariance Intersection%%%%
      xr_net_CI(:,t,i) = DCLCI(N, Gdyn(:,:,t), i, permute(xr_CI(:,t,:),[1 3 2]), permute(xr_IMU_CI(:,t,:),[1 3 2]), xTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      phir_net_CI(:,t,i) = DCLCI_phi(N, Gdyn(:,:,t), i, permute(phir_CI(:,t,:),[1 3 2]), permute(phir_IMU_CI(:,t,:),[1 3 2]), phiTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      
      
      %%%%%%%%%Covariance Union%%%%%%%%%%%
      xr_net_CU(:,t,i) = DCLCU(N, Gdyn(:,:,t), i, permute(xr_CU(:,t,:),[1 3 2]), permute(xr_IMU_CU(:,t,:),[1 3 2]), xTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      phir_net_CU(:,t,i) = DCLCU_phi(N, Gdyn(:,:,t), i, permute(phir_CU(:,t,:),[1 3 2]), permute(phir_IMU_CU(:,t,:),[1 3 2]), phiTg(:,t), 0.5*mu_mean(t,:), IMU_cov(t,:), mu_cov(t,:));
      
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:N
          xr_net_IMUf(:,t,i) = IMU_model(xr_IMUf(:,t,i),IMU_mean(t,i));  
          phir_net_IMUf(:,t,i) = IMU_model_phi(phir_IMUf(:,t,i),IMU_mean(t,i)); 
        end
   
        for i = 1:N
          xr_net_RSf(:,t,i) = xTg(:,t) - RangeSensing(xr_RSf(:,t,i),xTg(:,t),mu_mean(t,i));  
          phir_net_RSf(:,t,i) = phiTg(:,t) - RelPhiSens(phir_RSf(:,t,i),phiTg(:,t),mu_mean(t,i));
        end  
   
          xr_net(:,t,:) = xr_Est(:,t,:); 
          phir_net(:,t,:) = phir_Est(:,t,:); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
   % robot dynamics: pose sequence with DL-DCL    
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
   
   % robot dynamics: pose sequence with IMU only 
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
   
   
   % robot dynamics: pose sequence with RPSS only
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
   
   
   % robot dynamics: pose sequence with Kalman Fusion
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




   % robot dynamics: pose sequence with Covariance Intersection
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
   
   
   % robot dynamics: pose sequence with Covariance Union
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
   
end 
%% Calculation of the cumulative estimation losses for plots

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

%%%%%%%%%%%%%%%%%%%%%%%

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

% used when doing a Monte Carlo Sim.
for nrbt = 1:6

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

%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%
figure()
plot((1:T)*deltaT,CL_pos_nrbt,(1:T)*deltaT,CL_pos_nrbt_IMUf,(1:T)*deltaT,CL_pos_nrbt_RSf,(1:T)*deltaT,CL_pos_nrbt_KF,(1:T)*deltaT,CL_pos_nrbt_CI,(1:T)*deltaT,CL_pos_nrbt_CU)
set(gca,'FontSize',12)
xlabel('time (sec.)','FontSize',12,'Interpreter','Latex');
ylabel(sprintf('Cumulative Estimation Loss: $\\hat{L}_{t,%d}$',nrbt),'FontSize',12,'Interpreter','Latex');
title(sprintf('Cumul. Estim. loss of the ${%d}^{th}$ robot vs time',nrbt),'FontSize',12,'Interpreter','Latex');
legend('DL-DCL','No Comm., IMU only','No Comm., RPSS only','KF','CI','CU','Interpreter','Latex')
grid on
grid minor

end


figure()
plot((1:T)*deltaT,Tot_CL_pos./N,(1:T)*deltaT,Tot_CL_pos_IMUf./N,(1:T)*deltaT,Tot_CL_pos_RSf./N,(1:T)*deltaT,Tot_CL_pos_KF./N,(1:T)*deltaT,Tot_CL_pos_CI./N,(1:T)*deltaT,Tot_CL_pos_CU./N)
set(gca,'FontSize',12)
xlabel('time (sec.)','FontSize',12,'Interpreter','Latex');
ylabel('Total Cumulative Estimation Loss of all robots','FontSize',12,'Interpreter','Latex');
title("Total cumul. estim. loss of all robots vs time",'FontSize',12,'Interpreter','Latex');
legend('DL-DCL','No Comm., IMU only','No Comm., RPSS only','KF','CI','CU','Interpreter','Latex')
grid on
grid minor
%% Simulation Video

Tot_cumuL_fhat = Tot_CL_pos;
Tot_cumuL_fhatKF = Tot_CL_pos_KF;
Tot_cumuL_fhatCI = Tot_CL_pos_CI;
Tot_cumuL_fhatCU = Tot_CL_pos_CU;

%%% Simulation Video Creation %%%
n1 = 1;

myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)

figure('units','pixels','position',[0 0 1920 1080])

for t = 1:T
    
    if t >= 1 
       flag_clr = 1;    
    end
    if t >= ceil(T/6) 
       flag_clr = 2;    
    end
    if t >= ceil(2*T/6)
       flag_clr = 3;    
    end
    
    subplot(2,3,1)
    plot(xTg(1,1*n1:t),xTg(2,1*n1:t),'r','LineWidth', 0.001)
    hold on
    plot(xTg(1,t)+dT_safe(1)*cos(0:pi/180:2*pi),xTg(2,t)+dT_safe(1)*sin(0:pi/180:2*pi),'--r','LineWidth', 0.001)
    hold on
    plot(xr_KF(1,t,1),xr_KF(2,t,1),'x',xr_KF(1,t,2),xr_KF(2,t,2),'x',xr_KF(1,t,3),xr_KF(2,t,3),'x',xr_KF(1,t,4),xr_KF(2,t,4),'x',xr_KF(1,t,5),xr_KF(2,t,5),'x',xr_KF(1,t,6),xr_KF(2,t,6),'x',xTg(1,t),xTg(2,t),'s','LineWidth', 2, 'MarkerSize', 8)
    hold on 
    plot([xr_KF(1,t,1) xr_KF(1,t,1)+12*cos(phir_KF(:,t,1))],[xr_KF(2,t,1) xr_KF(2,t,1)+12*sin(phir_KF(:,t,1))],'--',[xr_KF(1,t,2) xr_KF(1,t,2)+12*cos(phir_KF(:,t,2))],[xr_KF(2,t,2) xr_KF(2,t,2)+12*sin(phir_KF(:,t,2))],'--',[xr_KF(1,t,3) xr_KF(1,t,3)+12*cos(phir_KF(:,t,3))],[xr_KF(2,t,3) xr_KF(2,t,3)+12*sin(phir_KF(:,t,3))],'--',[xr_KF(1,t,4) xr_KF(1,t,4)+12*cos(phir_KF(:,t,4))],[xr_KF(2,t,4) xr_KF(2,t,4)+12*sin(phir_KF(:,t,4))],'--',[xr_KF(1,t,5) xr_KF(1,t,5)+12*cos(phir_KF(:,t,5))],[xr_KF(2,t,5) xr_KF(2,t,5)+12*sin(phir_KF(:,t,5))],'--',[xr_KF(1,t,6) xr_KF(1,t,6)+12*cos(phir_KF(:,t,6))],[xr_KF(2,t,6) xr_KF(2,t,6)+12*sin(phir_KF(:,t,6))],'--','LineWidth', 0.1)
    hold off
    axis equal
    xlim([xTg(1,t)-20 xTg(1,t)+20])
    ylim([xTg(2,t)-20 xTg(2,t)+20])
    xlabel('$x_1$ (m)','FontSize',12,'Interpreter','Latex')
    ylabel('$x_2$ (m)','FontSize',12,'Interpreter','Latex')
    title('Kalman Fusion','FontSize',12,'Interpreter','Latex')

    subplot(2,3,2)
    plot(xTg(1,1*n1:t),xTg(2,1*n1:t),'r','LineWidth', 0.001)
    hold on
    plot(xTg(1,t)+dT_safe(1)*cos(0:pi/180:2*pi),xTg(2,t)+dT_safe(1)*sin(0:pi/180:2*pi),'--r','LineWidth', 0.001)
    hold on
    plot(xr_CI(1,t,1),xr_CI(2,t,1),'x',xr_CI(1,t,2),xr_CI(2,t,2),'x',xr_CI(1,t,3),xr_CI(2,t,3),'x',xr_CI(1,t,4),xr_CI(2,t,4),'x',xr_CI(1,t,5),xr_CI(2,t,5),'x',xr_CI(1,t,6),xr_CI(2,t,6),'x',xTg(1,t),xTg(2,t),'s','LineWidth', 2, 'MarkerSize', 8)
    hold on 
    plot([xr_CI(1,t,1) xr_CI(1,t,1)+12*cos(phir_CI(:,t,1))],[xr_CI(2,t,1) xr_CI(2,t,1)+12*sin(phir_CI(:,t,1))],'--',[xr_CI(1,t,2) xr_CI(1,t,2)+12*cos(phir_CI(:,t,2))],[xr_CI(2,t,2) xr_CI(2,t,2)+12*sin(phir_CI(:,t,2))],'--',[xr_CI(1,t,3) xr_CI(1,t,3)+12*cos(phir_CI(:,t,3))],[xr_CI(2,t,3) xr_CI(2,t,3)+12*sin(phir_CI(:,t,3))],'--',[xr_CI(1,t,4) xr_CI(1,t,4)+12*cos(phir_CI(:,t,4))],[xr_CI(2,t,4) xr_CI(2,t,4)+12*sin(phir_CI(:,t,4))],'--',[xr_CI(1,t,5) xr_CI(1,t,5)+12*cos(phir_CI(:,t,5))],[xr_CI(2,t,5) xr_CI(2,t,5)+12*sin(phir_CI(:,t,5))],'--',[xr_CI(1,t,6) xr_CI(1,t,6)+12*cos(phir_CI(:,t,6))],[xr_CI(2,t,6) xr_CI(2,t,6)+12*sin(phir_CI(:,t,6))],'--','LineWidth', 0.1)
    hold off
    axis equal
    xlim([xTg(1,t)-20 xTg(1,t)+20])
    ylim([xTg(2,t)-20 xTg(2,t)+20])
    xlabel('$x_1$ (m)','FontSize',12,'Interpreter','Latex')
    ylabel('$x_2$ (m)','FontSize',12,'Interpreter','Latex')
    title('Covariance Intersection','FontSize',12,'Interpreter','Latex')
  
    subplot(2,3,3)
    plot(deltaT*(1:t),Tot_cumuL_fhat(1:t)./N,deltaT*(1:t),Tot_cumuL_fhatKF(1:t)./N,deltaT*(1:t),Tot_cumuL_fhatCI(1:t)./N,deltaT*(1:t),Tot_cumuL_fhatCU(1:t)./N,'LineWidth',1.2)
    set(gca,'FontSize',12)
    xlabel('time (sec.)','FontSize',12,'Interpreter','Latex');
    ylabel('Avg. Cumul. Loss','FontSize',12,'Interpreter','Latex');
    title("Avg. cumul. loss per robot vs time",'FontSize',12,'Interpreter','Latex');
    legend('DL-DCL','KF','CI','CU','Interpreter','Latex')
    xlim([0 deltaT*T])
    ylim([0 90])
    grid on
    grid minor
    
    subplot(2,3,4)
    plot(xTg(1,1*n1:t),xTg(2,1*n1:t),'r','LineWidth', 0.001)
    hold on
    plot(xTg(1,t)+dT_safe(1)*cos(0:pi/180:2*pi),xTg(2,t)+dT_safe(1)*sin(0:pi/180:2*pi),'--r','LineWidth', 0.001)
    hold on
    plot(xr(1,t,1),xr(2,t,1),'x',xr(1,t,2),xr(2,t,2),'x',xr(1,t,3),xr(2,t,3),'x',xr(1,t,4),xr(2,t,4),'x',xr(1,t,5),xr(2,t,5),'x',xr(1,t,6),xr(2,t,6),'x',xTg(1,t),xTg(2,t),'s','LineWidth', 2, 'MarkerSize', 8)
    hold on 
    plot([xr(1,t,1) xr(1,t,1)+12*cos(phir(:,t,1))],[xr(2,t,1) xr(2,t,1)+12*sin(phir(:,t,1))],'--',[xr(1,t,2) xr(1,t,2)+12*cos(phir(:,t,2))],[xr(2,t,2) xr(2,t,2)+12*sin(phir(:,t,2))],'--',[xr(1,t,3) xr(1,t,3)+12*cos(phir(:,t,3))],[xr(2,t,3) xr(2,t,3)+12*sin(phir(:,t,3))],'--',[xr(1,t,4) xr(1,t,4)+12*cos(phir(:,t,4))],[xr(2,t,4) xr(2,t,4)+12*sin(phir(:,t,4))],'--',[xr(1,t,5) xr(1,t,5)+12*cos(phir(:,t,5))],[xr(2,t,5) xr(2,t,5)+12*sin(phir(:,t,5))],'--',[xr(1,t,6) xr(1,t,6)+12*cos(phir(:,t,6))],[xr(2,t,6) xr(2,t,6)+12*sin(phir(:,t,6))],'--','LineWidth', 0.1)
    hold off
    axis equal
    xlim([xTg(1,t)-20 xTg(1,t)+20])
    ylim([xTg(2,t)-20 xTg(2,t)+20])
    xlabel('$x_1$ (m)','FontSize',12,'Interpreter','Latex')
    ylabel('$x_2$ (m)','FontSize',12,'Interpreter','Latex')
    title('DL-DCL','FontSize',12,'Interpreter','Latex')
    
    subplot(2,3,5)
    plot(deltaT*(1:t),mu_mean(1:t,1).*(1/4),'o',deltaT*(1:t),mu_mean(1:t,2).*(2/4),'o',deltaT*(1:t),mu_mean(1:t,3).*(3/4),'o',deltaT*(1:t),mu_mean(1:t,4).*(4/4),'o',deltaT*(1:t),mu_mean(1:t,5).*(5/4),'o',deltaT*(1:t),mu_mean(1:t,6).*(6/4),'o')
    xlim([0 deltaT*T])
    ylim([0.8 6.2])
    set(gca,'FontSize',12)
    xlabel('time (sec.)','FontSize',12,'Interpreter','Latex')
    ylabel('$i^{th}$ robot','FontSize',12,'Interpreter','Latex')
    title('RPSS temporary failure vs time','FontSize',12,'Interpreter','Latex')
    %legend('Rbt1','Rbt2','Rbt3','Rbt4','Rbt5','Rbt6','Interpreter','Latex')
    grid on 
    
    subplot(2,3,6)
    hg = plot(graph(Gdyn(:,:,t)-diag(ones(1,6))));
    hg.XData = [1 1+sqrt(2) sqrt(2)+2 sqrt(2)+1 1 0];
    hg.YData = [0 0 1 1+sqrt(2) 1+sqrt(2) sqrt(2)];
    hg.LineWidth = 2;
    hg.NodeColor = [0 0.4470 0.7410;
                    0.8500 0.3250 0.0980;
                    0.9290 0.6940 0.1250;
                    0.4940 0.1840 0.5560;
                    0.4660 0.6740 0.1880;
                    0.3010 0.7450 0.9330];
                if flag_clr == 1
                  hg.NodeColor(IMU_ord(1),:) = [0 0 0];
                end
                if flag_clr == 2
                  hg.NodeColor(IMU_ord(1),:) = [0 0 0];
                  hg.NodeColor(IMU_ord(2),:) = [0 0 0];
                end
                if flag_clr == 3
                  hg.NodeColor(IMU_ord(1),:) = [0 0 0];
                  hg.NodeColor(IMU_ord(2),:) = [0 0 0];
                  hg.NodeColor(IMU_ord(3),:) = [0 0 0];
                end          
    hg.MarkerSize = 6;
    title('Communication/Interaction Network')
    xlabel('Nodes are robots (black node = failed IMU/FLT)')
    ylabel('Edges are communication links')
    pause(0.1)
    
    sgtitle(['Time = ',num2str(t*deltaT,'%4.1f'),' sec.']) 
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)