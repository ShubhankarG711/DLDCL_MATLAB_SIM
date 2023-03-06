function [phir_IMU] = IMU_model_phi(phir_curr,IMU_mean)
% add bias to the actual heading angle
phir_IMU = wrapToPi(phir_curr + ((10*pi/180)/(0.9)).*IMU_mean.*0.5);
end
