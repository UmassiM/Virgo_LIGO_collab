function [hfdf,reduced_gridk] = LAB_hough_and_peak_noantsour(peaks,Tfft,freq_band,this_n,kn)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

peaks(1,:)=peaks(1,:)/86400;
ref_time=0;
[gridk,~]=LAB_long_transient_grid_k(Tfft,freq_band,[1/(2*Tfft)^2 1/Tfft^2],this_n);

hm_job.minf=freq_band(1);% minimum frequency to do the Hough on
hm_job.maxf=freq_band(2);% maximum frequency to do the Hough on, the value of this
hm_job.df=1/Tfft; %step in frequency
hm_job.frenh=1;

kmin=kn/1.15; % max/min values you want to consider
kmax=kn*1.15;
[~, ind1]=min(abs(kmin-gridk));
[~, ind2]=min(abs(kmax-gridk));
if length(gridk)>800    
    reduced_gridk=gridk(ind1:ind2);    
else
    reduced_gridk=gridk; % zooming in on portion of grid in k that is the near the source
end

[hfdf,grid_x]=GENERALIZED_transients_hough_nogd(peaks,hm_job,this_n,ref_time,reduced_gridk);

end

