function [dx_host_obsvd] = RangeSensing(xr_host,xr_obsvd,mu_mean)
% add bias to the RPSS relative position estimate
dx_host_obsvd = (xr_obsvd - xr_host) + mu_mean.*0.5;
end

