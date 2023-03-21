function [Y,NN1,NN2]=thorus_recenter(y,inv)
% recenter a thoroidal image (for fourier transforms)
%
%   y     matrix or image structure
%   inv   inverse operation

% Sapienza Università di Roma
% Laboratorio di Segnali e Sistemi II
% Author: Sergio Frasca - 2018

if ~exist('inv','var')
    inv=0;
end

if isstruct(y)
    [N1,N2]=size(y.Z);
else
    [N1,N2]=size(y);
end

NN1=ceil(N1/2);
NN2=ceil(N2/2);
    
if inv > 0
    Y=rota2(y,-NN1,-NN2);
else
    Y=rota2(y,NN1,NN2);
end