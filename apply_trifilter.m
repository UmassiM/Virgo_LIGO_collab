function out=apply_trifilter(mapin,m1,m2,cut)
%
%
%
% %
%
out=mapin;
if isstruct(mapin)
    map0=mapin.Z;
else
    map0=mapin;
end
[nx,ny]=size(map0);
M1=-ny/(nx*m1);
M2=-ny/(nx*m2);

f0=generate_trifilter(nx,ny,M1,M2,cut);
% NN1=floor(nx/2);
% NN2=floor(ny/2);
% f1=f0;
% f0(1:NN1,:)=f1(NN1+1:nx,:);
% f0(NN1+1:nx,:)=f1(1:NN1,:);
% f1=f0;
% f0(:,1:NN2)=f1(:,NN2+1:ny);
% f0(:,NN2+1:ny)=f1(:,1:NN2);
[f0,NN1,NN2]=thorus_recenter(f0);
mapfft=fft2(map0);
mapfiltered=ifft2(mapfft.*f0);
if isstruct(mapin)
    out.Z=real(mapfiltered);
else
    out=real(mapfiltered);
end

