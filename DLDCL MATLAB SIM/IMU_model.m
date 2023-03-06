function [xr_IMU] = IMU_model(xr_curr,IMU_mean)
% add bias to the actual position
xr_IMU = xr_curr + IMU_mean*0.5;
end

