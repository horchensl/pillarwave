%make taper
%
% usage tap=mktaper(l,n);
%
% input l : length of array
%       n : # taper points;
%
%
%
%
%
%--------------------------------------------

function tap=mktaper(l,n);

if (2*n) <= l
   x=0:pi/(n-1):pi; 
   x=(1+cos(x))/2;
   m=n:-1:1;
   tap=[x(m),ones(1,l-2*n),x];
else
   tap=[];
end

