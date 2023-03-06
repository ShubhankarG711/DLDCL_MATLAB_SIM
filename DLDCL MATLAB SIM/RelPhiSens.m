function [dphi_host_obsvd] = RelPhiSens(phir_host,phir_obsvd,mu_mean)
% add bias to the RPSS relative heading angle estimate
dphi_host_obsvd = wrapToPi((phir_obsvd - phir_host) + ((10*pi/180)/(0.9)).*mu_mean.*0.5);
end
