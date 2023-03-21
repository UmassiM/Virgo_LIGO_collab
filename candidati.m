function [x_candidato,k_candidato] = candidati(x0,kn,gridx,gridk)
[row,col]=min(abs(kn-gridk(:)));
[row2,col2]=min(abs(x0-gridx(:)));
x_candidato=gridx(col2);
k_candidato=gridk(col);


end

