function image_plot(im_str,scale,op)
% image_plot  plot image 
%
%    im_str   input image structure
%    scale    [min max] of z scale or 0
%    op       0 normal, 1 log, 2 sqrt
%
%  im_str.Z  (nx,ny)
%  im_str.inix
%  im_str.dx
%  im_str.iniy
%  im_str.dy

% Sapienza Università di Roma
% Laboratorio di Segnali e Sistemi II
% Author: Sergio Frasca - 2018

if isstruct(im_str)
    M=im_str.Z;
    [n1,n2]=size(M);
    x=im_str.inix+(0:n1-1)*im_str.dx;
    y=im_str.iniy+(0:n2-1)*im_str.dy;
else
    M=im_str;
    [n1,n2]=size(M);
    x=1:n1;
    y=1:n2;
end

s2.M=M;
s2.x=x;
s2.y=y;

if exist('op','var')
    switch op
        case 1
            M=log10(abs(M));
        case 2
            M=sqrt(abs(M));
    end
else
    op=0;
end

if ~exist('scale','var') | scale == 0
    M1=M(:);
    i=find(M1 ~= 0);
    M1=M1(i);
    scale=[min(M1) max(M1)];
end

h0scra = figure('Color',[1 1 0.6], ...
   'Name','gd2 Image');

if op == 0
    him=imagesc(x,y,M',scale);
else
    him=imagesc(x,y,M');
end

colorbar
set(gca,'YDir','normal')

set(gcf,'UserData',s2)
strget='s2=get(gcf,''UserData'');M=s2.M;x=s2.x;y=s2.y;';

if op == 0
    dismenuh=uimenu(h0scra,'label','Display','UserData',him);

    uimenu(dismenuh,'label','Normal','callback',...
       [strget 'imagesc(x,y,M'',scale);grid on;zoom on;colorbar,set(gca,''YDir'',''normal'')']);
    uimenu(dismenuh,'label','Log','callback',...
       [strget 'imagesc(x,y,log10(M''));grid on;zoom on;colorbar,set(gca,''YDir'',''normal'')']);
    uimenu(dismenuh,'label','Square Root','callback',...
       [strget 'imagesc(x,y,sqrt(abs(M'')));grid on;zoom on;colorbar,set(gca,''YDir'',''normal'')']);
end

opmenuh=uimenu(h0scra,'label','Operations');

uimenu(opmenuh,'label','x projection','callback',...
   [strget  'mx_gd2m=mean(M,2);figure;plot(x,mx_gd2m),grid on,title(''x projection''),xlim([min(x) max(x)])' ...
   ';hold on,if length(x) < 256; plot(x,mx_gd2m,''r.''); end']);
uimenu(opmenuh,'label','y projection','callback',...
   [strget 'my_gd2m=mean(M,1);figure;plot(y,my_gd2m),grid on,title(''y projection''),xlim([min(y) max(y)])' ...
   ';hold on,if length(y) < 256; plot(y,my_gd2m,''r.''); end']);
uimenu(opmenuh,'label','x projection (log)','callback',...
   [strget 'mx_gd2m=mean(M,2);figure;semilogy(x,mx_gd2m),grid on,title(''x projection''),xlim([min(x) max(x)])']);
uimenu(opmenuh,'label','y projection (log)','callback',...
   [strget 'my_gd2m=mean(M,1);figure;semilogy(y,my_gd2m),grid on,title(''y projection''),xlim([min(y) max(y)])']);

uimenu(opmenuh,'label',' ');
uimenu(opmenuh,'label','Selection','callback',...
   [strget 'p=ginput(2),[kx ky]=coord2pix_gd2II(x,y,p(:,1),p(:,2));'...
   'kx=[min(kx) max(kx)],ky=[min(ky) max(ky)],'...
   'g3.mat=M(kx(1):kx(2),ky(1):ky(2));g3.x=x(kx(1):kx(2));g3.y=y(ky(1):ky(2));M=g3.mat;x=g3.x;y=g3.y;g3,image_gd2(g3)']);
