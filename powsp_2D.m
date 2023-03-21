function [outsp,outft]=powsp_2D(in,res,rec,icplot)
% 2-D power spectrum (for 2D image structure)
%
%   in      input 2D structure 
%   res     [resx resy] resolution enhancement; < 0 no mean subtraction
%   rec     > 0 -> recenter the spectrum
%   icplot  > 0 plot; def 0

% Sapienza Università di Roma
% Laboratorio di Segnali e Sistemi II
% Author: Sergio Frasca - 2018

if ~exist('icplot','var')
    icplot=0;
end
nomean=1;
if res(1) < 0
    res=-res;
    nomean=0;
end
    
if length(res) == 1
    res=[res res];
end

if ~exist('rec','var')
    rec=0;
end

y=in.Z;
if nomean > 0
    y=y-mean(y(:));
end
rec=abs(rec);

dx1=in.dx;
dx2=in.dy;

[n1,n2]=size(y);

% N1=round(n1*res(1)/2)*2;
% N2=round(n2*res(2)/2)*2;
N1=round(n1*res(1));
N2=round(n2*res(2));

Y=zeros(N1,N2);

Y(1:n1,1:n2)=y;

Y=fft2(Y);
if rec == 0
    outft=Y;
else
    outft=0;
end
Y=(abs(Y).^2)/(N1*N2);

df1=1/(dx1*N1);
df2=1/(dx2*N2);
ini1=0;
ini2=0;

if rec > 0
%     NN1=ceil(N1/2);
%     NN2=ceil(N2/2);
%     y=Y;
%     Y(1:NN1,:)=y(NN1+1:N1,:);
%     Y(NN1+1:N1,:)=y(1:NN1,:);
%     y=Y;
%     Y(:,1:NN2)=y(:,NN2+1:N2);
%     Y(:,NN2+1:N2)=y(:,1:NN2);
    [Y,NN1,NN2]=thorus_recenter(Y);
    ini1=-NN1*df1;
    ini2=-NN2*df2;
end

outsp.Z=Y;
outsp.dx=df1;
outsp.dy=df2;
outsp.inix=ini1;
outsp.iniy=ini2; 

if icplot > 0
    image_plot(outsp,0);
end

% out=edit_gd2(out,'dx',df1,'dx2',df2,'ini',ini1,'ini2',ini2);