function [] = critical(hfdf,gridk, gridx, a)
%ciò che sta qui è sicuramente spiegato in "sectors"
hfdf0=hfdf(1:floor(2464/4.77),:);
x=size(hfdf0);
x0=x(1);
y0=x(2);
x=floor(x0/5);
y=floor(y0/5);
figure,uimagesc(gridx(1:floor(2464/4.77)),gridk,hfdf0'),colorbar,axis xy, title(compose('quarto'))
i=3; %queste sono le coordinate del quadrante contenente il max nel caso del generated long transient con mappa divisa per 4.77, andrebbero trovate a mano
j=3;

pixels=hfdf0(1+(i-1)*x:i*x, 1+(j-1)*y:j*y);
[row,col]=find(pixels==max(pixels(:)));
pos=[(i-1)*x+row(1),(j-1)*y+col(1)];
row=pos(1);
col=pos(2);
cclose=hfdf0(max(1,row-floor(a*3)):min(x0,row+floor(a*3)),max(1,col-floor(a*1.5)):min(y0,col+floor(a*1.5)));
figure,uimagesc(gridx(1+(i-1)*x:i*x),gridk(1+(j-1)*y:j*y),pixels'),colorbar,axis xy, title(compose('quadrante'))
figure,uimagesc(gridx(max(1,row-floor(a*3)):min(x0,row+floor(a*3))),gridk(max(1,col-floor(a*1.5)):min(y0,col+floor(a*1.5))),cclose'),colorbar,axis xy, title(compose('intorno'))

size(pixels)
size(cclose)

%trovo il valore critico sia per il quadrante, sia per cclose
m=max(pixels(:));
avg=mean(pixels(:));
dev=std(pixels(:));
crit=(m-avg)/dev

m=max(cclose(:));
avg=mean(cclose(:));
dev=std(cclose(:));
crit=(m-avg)/dev

end