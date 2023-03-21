function plot_triplets(x,y,z,marker,colorm)
%PLOT_TRIPLETS plots many points with different colors
%              faster than plot_manypoints
%
%          plot_triplets(x,y,z,marker,colorm)
%
%      x,y,z      triplets; x,y coordinates, z amplitude
%      marker     'x', 'o','+','.','<','>','^','v','s','d',...
%      colormap   default 'cool'

% Version 2.0 - August 2006
% Part of Snag toolbox - Signal and Noise for Gravitational Antennas
% by Sergio Frasca - sergio.frasca@roma1.infn.it
% Department of Physics - Universita` "La Sapienza" - Rome

fig = get(groot,'CurrentFigure');
if ~isempty(fig)
    figure;
end
mi=min(z);
ma=max(z);

mami=ma-mi;
if mami <= 0
    mami=1
    disp(' Attention ! mami = 1')
end

if ~exist('marker','var')
    marker='.';
end

if ~exist('colorm','var')
    colorm='jet';
end

cm=colormap(colorm); 
% figure
nc=length(cm);
zz=floor(nc*0.9999*(z-mi)/mami+1);

for i = 1:nc
    col=cm(i,:);
    plot(x(zz==i),y(zz==i),marker,'MarkerSize',4,'MarkerFaceColor',col,'MarkerEdgeColor',col),hold on
end

caxis([min(z) max(z)])
colormap(colorm),colorbar
grid on