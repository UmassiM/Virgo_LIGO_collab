function y=ps_window(typ,len,par)
% computes some windows for power spectrum estimates
%
%         y=ps_window(typ,len,par)
%
%  type    a number or string as 
%            'bartlett'  1 Bartlett window
%            'hanning'   2 Hanning window
%            'flatcos'   3 flat-top, cosine edge
%            'tukey'     4 Tukey window 
%            'gauss'     5 gaussian (3 sigma)
%            'no'        0 no window (flat)
%  len     length
%  par     parameter(s) (if needed)

% Sapienza Università di Roma
% Laboratorio di Segnali e Sistemi II
% Author: Sergio Frasca - 2018

y=ones(1,len);
len2=len/2;
len4=ceil(len/4);

if isnumeric(typ)
    switch typ
        case 1
            typ='bartlett';
        case 2
            typ='hanning';
        case 3
            typ='flatcos';
        case 4
            typ='tukey';
        case 5
            typ='gauss';
    end
end

switch typ
case 'bartlett'
   y(1:len2)=(1:len2)./len2;
   y(len:-1:len2+1)=y(1:len2);
   y=y*sqrt(3);
case 'hanning'
   y(1:len2)=(1-cos((1:len2)*pi/len2));
   y(len:-1:len2+1)=y(1:len2);
   y=y*sqrt(2/3);
case 'flatcos'
   y(1:len4)=(1-cos((1:len4)*pi/len4))/2;
   y(len:-1:len-len4+1)=y(1:len4);
   y=y*sqrt(len/sum(y.^2));
case 'tukey'  % par = alpha
   lenx=round(par*len/2);
   y(1:lenx)=(1-cos((1:lenx)*pi/lenx))/2;
   y(len:-1:len-lenx+1)=y(1:lenx);
   y=y*sqrt(len/sum(y.^2));
case 'gauss'
   dx=6/(len-1);
   x=-3:dx:3;
   y=exp(-x.^2/2)/sqrt(2*pi);
   y=y*sqrt(len/sum(y.^2));
end
