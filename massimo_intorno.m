function out = massimo_intorno(hfdf,x0,kn,gridx,gridk)
[row,col]=min(abs(kn-gridk(:)));
[row2,col2]=min(abs(x0-gridx(:)));
out=max(hfdf(col2-2:col2+2,col-2:col+2),[],'all');

end

