function [find_k,find_x]=difference_gridx_gridk(kn,gridk,x0,gridx)
find_k=min(abs(kn-gridk(:)));
find_x=min(abs(x0-gridx(:)));
% [row,col]=min(abs(kn-gridk(:)))
% [row2,col2]=min(abs(x0-gridx(:)))
end




