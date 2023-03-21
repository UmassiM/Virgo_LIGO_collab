function image_gd2_nonunigrid(hfdf,y,n,x)
% image_gd2_nonunigrid  image of a gd2 that has at least 1 non uniform grid
% 
%
% x and y can be nonuniformly spaced !!
%
%    g2       input gd2
%    scale    [min max] of z scale or 0
%    op       0 normal, 1 log, 2 sqrt; op=0 for me

% Version 2.0 - November 2010
% Part of Snag toolbox - Signal and Noise for Gravitational Antennas
% by Sergio Frasca - sergio.frasca@roma1.infn.it
% Department of Physics - Universit� "Sapienza" - Rome
% 
% if isstruct(g2)
%     M=g2.mat;
%     x=g2.x;
%     y=g2.y;
% elseif isreal(g2)
%     [i1 i2]=size(g2);
%     M=g2;
%     x=1:i1;
%     y=1:i2;
% else
%     M=y_gd2(g2);
%     x=x_gd2(g2);
%     y=x2_gd2(g2);
% end

if exist('n','var')
    if n~=1
        outt=x_gd2(hfdf).^(1/(1-n));
        hfdf=edit_gd2(hfdf,'x',sort(outt),'y',flip(y_gd2(hfdf)));
    else
        outt=exp(x_gd2(hfdf));
        hfdf=edit_gd2(hfdf,'x',outt,'y',flip(y_gd2(hfdf)));
    end
end


M=y_gd2(hfdf);
if ~exist('x','var') %basically the case where grid in x is uniform 
    x=x_gd2(hfdf);
end

if ~exist('x','var') & ~exist('y','var')
    image_gd2(hfdf);
    return;
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
    him=uimagesc(x,y,M',scale);
else
    him=uimagesc(x,y,M');
end
colormap(jet(256))%%%%%how to adjust colorbar
colorbar
set(gca,'YDir','normal')
set(gca, 'FontSize', 14)
set(gcf,'UserData',s2)
strget='s2=get(gcf,''UserData'');M=s2.M;x=s2.x;y=s2.y;';

if op == 0
    dismenuh=uimenu(h0scra,'label','Display','UserData',him);

    uimenu(dismenuh,'label','Normal','callback',...
       [strget 'uimagesc(x,y,M'',scale);grid on;zoom on;colorbar,set(gca,''YDir'',''normal'')']);
    uimenu(dismenuh,'label','Log','callback',...
       [strget 'uimagesc(x,y,log10(M''));grid on;zoom on;colorbar,set(gca,''YDir'',''normal'')']);
    uimenu(dismenuh,'label','Square Root','callback',...
       [strget 'uimagesc(x,y,sqrt(abs(M'')));grid on;zoom on;colorbar,set(gca,''YDir'',''normal'')']);
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
