function  out= Pixel(hfdf,x0,kn,gridx,gridk)
[row,col]=min(abs(kn-gridk(:)));
[row2,col2]=min(abs(x0-gridx(:)));
 
out=hfdf(max(col2-100,1):col2+100, col-30:col+30);%costruisco la sottomappa

% image_plot(out);

end

