function [ID_nrst,dist_nrst] = NrstNeighDist(rbt_ID, xr_all)
% find out the robot which is spatially nearest to robot with ID as rbt_ID   

xr_curr = xr_all(:,rbt_ID);
dist_all = sqrt(sum((xr_curr - xr_all).^2,1));

[dist_srt,ID_srt] = sort(dist_all);

ID_nrst = ID_srt(2);
dist_nrst = dist_srt(2);
end