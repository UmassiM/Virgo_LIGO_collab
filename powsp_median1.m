function [pow_med,f_med]=powsp_median(one_side_fft_abs,f,bandw,nover)
% 
% 
% 
% 
pow_data=abs(one_side_fft_abs).^2;
medlength=bandw;
medover=floor(medlength/nover);
N=length(pow_data);
Nmeds=floor(N/medover);
pow_med=zeros(1,Nmeds);
f_med=zeros(1,Nmeds);
i_f=1;
while (i_f+1)*medover<N
    pow_med(i_f)=median(pow_data(((i_f-1)*medover+1):((i_f+1)*medover)),'omitnan');
    f_med(i_f)=f(i_f*medover);
    i_f=i_f+1;
end
pow_med=sqrt(pow_med/log(2));