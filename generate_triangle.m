function out=generate_triangle(nx,ny,m1,m2,ycut)
%
%
%
%
%
map=zeros(nx,ny);
M1=min(m1,m2);
M2=max(m1,m2);
if ycut>0 && ycut<1
    ytop=floor(ny/2)+1+round(ycut*ny/2);
else
    ytop=ny;
end
if rem(nx,2)==0
    i01=floor(nx/2);
    i02=i01+1;
%     i01=floor(nx/2)+1;
%     i02=i01;
else
    i01=floor(nx/2)+1;
    i02=i01;
%     i01=floor(nx/2);
%     i02=i01;
end
% if rem(ny,2)==0
%     j01=floor(ny/2)+1;
%     j02=j01;
% else
%     j01=floor(ny/2);
%     j02=j01;
% end
j01=floor(ny/2)+1;
j02=j01;
for j=j01:ytop
    if abs(M1)==Inf
        i1=i01;
    else
        i1=(j-j01)/M1+i01;
    end
    if abs(M2)==Inf
        i2=i02;
    else
        i2=(j-j02)/M2+i02;
    end
    i11=round(min(i1,i2));
    i22=round(max(i1,i2));
    imin=max(1,i11);
    imax=min(i22,nx);
    map(imin:imax,j)=1;
end
mirror=zeros(nx,j01-1);
if rem(ny,2)==0
    mirror(:,:)=map(nx:-1:1,j01:ny);
else
    mirror(:,:)=map(nx:-1:1,(j01+1):ny);
end
map(:,1:(j01-1))=mirror(:,(j01-1):-1:1);
out=map;
