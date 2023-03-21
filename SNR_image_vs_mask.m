function out=SNR_image_vs_mask(map,mask)
%
%
%
values_in=map(mask>=1);
values_out=map(mask<1);
out=(mean(values_in)-mean(values_out))/std(values_out);
