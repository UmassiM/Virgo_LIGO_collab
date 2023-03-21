function out = Intorno(x0,kn,gridx,gridk,hfdf)
[row,col]=min(abs(kn-gridk(:)));
[row2,col2]=min(abs(x0-gridx(:)));
out=hfdf(col2-2:col2+2,col-2:col+2);

end

