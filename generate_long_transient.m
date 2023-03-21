function [amps,ff]=generate_long_transient(f0,kn,h0,n,tt)
%
%
%
%
% pc=3.086*10^16;
tau=1/(f0^(n-1)*kn*(n-1));

ff=f0.*((1+tt./tau).^(1/(1-n)));
% df=-(f0^(1-n))*(ff.^n)/((n-1)*tau);
% phase=2*pi*mod(tau*f0*((1-n)/(2-n)).*((1+tt./tau).^((2-n)/(1-n))-1),1);
% h0=1.8*10^(-24)*20*pc*(10^6).*((ff./1000).^3)./(DL*pc*(10^6));
if n<6.8 && n~=11/3
    amps=h0*(ff.^2)/f0^2;
elseif n>=6.8 && n<7.2
    amps=h0*(ff.^3)/f0^3;
elseif n==11/3
    amps=h0*ff.^(2/3)/f0^(2/3);
end
% hsig=amps.*sin(phase);