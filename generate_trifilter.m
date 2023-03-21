function out=generate_trifilter(nx,ny,m1,m2,cut)
%
%
%   cut: relative frequency cut (0<>1)
%         if 2-comp, [cutx cuty]
%
map=zeros(nx,ny);
M1=min(m1,m2);
M2=max(m1,m2);


i01=ceil(nx/2)+1;
i02=i01;
j01=ceil(ny/2)+1;
j02=j01;

for j=1:ny
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

if length(cut)==1
    cut(2)=cut;
end
if cut(1)>0 && cut(1)<1
    cutx=round((1-cut(1))*nx/2);
    map(1:cutx,:)=0;
    map((nx-cutx+1):nx,:)=0;
end
if cut(2)>0 && cut(2)<1
    cuty=round((1-cut(2))*ny/2);
    map(:,1:cuty)=0;
    map(:,(ny-cuty+1):ny)=0;
end

out=map;
