function [] = sectors(l, hfdf)

x=size(hfdf);
hfdf1=hfdf(:,1:floor(x(2)/4.77)); %divisione opzionale, è comoda per scartare le alte frequenze
x=size(hfdf1);
x0=x(1);
y0=x(2);
x=floor(x0/5);
y=floor(y0/5);

for k=[1:1:l] %l serve ad iterare vari valori della dimensione di close tra a/l ed a
    a=floor(y/2/l*k); %il valore massimo (k=l) fa in modo di prendere tutto il quadrante (o quasi, col floor)
    crit=zeros(5,5);
    for i=[1:1:5] %cicli per esaminare ogni quadrante
        for j=[1:1:5]
            pixels=hfdf1(1+(i-1)*x:i*x, 1+(j-1)*y:j*y); %quadrante di interesse
            [row,col]=find(pixels==max(pixels(:))); %trovo il max nel quadrante (in caso di più coordinate si usa la prima)
            pos=[(i-1)*x+row(1),(j-1)*y+col(1)]; %trovo le coordinate nel riferimento completo
            row=pos(1); %le sostituisco a quelle del rif del quadrante
            col=pos(2);
            %ho detto abbastanza di close nella mail: sono parametri
            %trovati su misura della mappa, ma si possono sostituire con
            %costanti di proporzionalità in funzione di x ed y
            cclose=hfdf1(max(1,row-floor(a*3)):min(x0,row+floor(a*3)),max(1,col-floor(a*1.5)):min(y0,col+floor(a*1.5)));
            m=max(cclose(:));
            avg=mean(cclose(:));
            %med=median(close(:));
            dev=std(cclose(:));
            crit(i,j)=(m-avg)/dev;
        end
    end
    crit
end

end