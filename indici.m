function [row,col,row2,col2] = indici(x0,kn,gridx,gridk)

[row,col]=min(abs(kn-gridk(:)));
[row2,col2]=min(abs(x0-gridx(:)));

end

