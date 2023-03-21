function [pow_med]=powsp_median1(pow_data,bandw,nover)
% 
% 
% 
% 
medlength=bandw;
medover=floor(medlength/nover);
N=length(pow_data);
Nmeds=floor(N/medover);
pow_med=zeros(1,Nmeds);
i_f=1;


while (i_f+1)*medover<N
    pow_med(i_f)=median(pow_data(((i_f-1)*medover+1):((i_f+1)*medover)),'omitnan');
    i_f=i_f+1;  
end

pow_med=pow_med/log(2);

