function out=rota2(in,k1,k2)
% rotation for images and matrices in the 2 dimensions (on thorus)
%
%     out=rota2(in,k1,k2)
%
%   in      image structure or matrix
%   k1,k2   rotation bins

% Sapienza Università di Roma
% Laboratorio di Segnali e Sistemi II
% Author: Sergio Frasca - 2018

icstr=0;
if isstruct(in)
    out=in;
    in=in.Z;
    icstr=1;
end

[N1,N2]=size(in); 
y=in;

if k1 < 0
    k1=N1+k1;
end
if k2 < 0
    k2=N2+k2;
end
Y(1:N1-k1,:)=y(k1+1:N1,:);
Y(N1-k1+1:N1,:)=y(1:k1,:);
y=Y;
Y(:,1:N2-k2)=y(:,k2+1:N2);
Y(:,N2-k2+1:N2)=y(:,1:k2);

if icstr
    out.Z=Y;
else
    out=Y;
end