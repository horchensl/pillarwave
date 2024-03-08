function z = cor(x,y,ipad,ifull)
%COR	Correlation via frequency domain.
%	COR(X,Y) returns the cross-correlation product of X and Y
%TSMAT	via multiplication in the complex (=frequency) domain.
%TLBX	If X and\or Y are real, FFT's are done.
%	If both X and Y are real the result is also returned
%	as real.
%	For real vectors of unequal lengths, the shortest is padded
%	with zeros prior to the FFT, and the result is returned as
%	long as the longest.
%
%	Correlation via the frequency domain yields the so-called
%	circular correlation which, when padded with extra zeros,
%	is identical upto machine accuracy to the time domain
%	correlation as calculated with XCORR but for longer signals
%	very much faster.
%	The circular correlation normally returns only the positive
%	lags, i.e. from 0 to N-1 for length N signals.
%	For padding with N zeroes during calculation, where N is
%	the length of the longest signal, use COR(X,Y,IPAD) with
%	IPAD = 1.
%	For a full correlation including the negative time lags,
%	i.e. running from -N+1 ... 0 ... N-1, use COR(X,Y,IPAD,IFULL)
%	with IPAD and IFULL both 1. This result should be equal to
%	XCORR(X,Y).
%	Also note that XCORR in Matlab 5.x is erroneously calculated
%	time reversed. This has been corrected in Matlab 6.x.
%
%	See e.g. "Bendat, J.S., and A.G. Piersol. Random Data: Analysis
%	and Measurement Procedures. New York: John Wiley & Sons, 1971."
%	for algorithm ideas and correlation function properties.
%
%	See also XCORR (Signal Processing Toolbox).

%	Uilke Stelwagen, October 1991, November 2000.
%	Copyright (C) 1991, 2000 Institute of Applied Physics, TNO-TPD,
%	The Netherlands.

if nargin==1
   y = x; ipad = 0; ifull = 0;
elseif nargin==2
   ipad = 0; ifull = 0;
elseif nargin==3
   ifull = 0;
elseif nargin~=4
   help cor
   return
end

ixb=0; iyb=0;
nx=size(x,1); ny=size(y,1);
ns = max([nx ny]);
if ipad
   nfft = 2*ns;
else
   nfft = ns;
end
if(isreal(x))
  x = fft(x,nfft);
  ixb=1;
end
if(isreal(y))
  y = fft(y,nfft);
  iyb=1;
end

z = x.*conj(y);     % yields positive time lags
if ifull
   zn = y.*conj(x); % yields negative time lags because
end                 % Rxy(-tau) = Ryx(tau)*, e.g. see Bendat and Piersol, pp. 29.
   
if(ixb==1 & iyb==1)
   z = real(ifft(z));
   if ifull
      zn = real(ifft(zn));
   end
end

if ipad
   z = z(1:ns);
   if ifull
      zn = zn(1:ns);
   end
end

if ifull
   [m,n] = size(z);
   if m>n
      z = [(zn(ns:-1:2)); z]; % column vector
   else
      z = [(zn(ns:-1:2))  z]; % row vector
   end
end

if(ixb~=iyb)
  disp(' Warning: input vectors in different domains')
  disp(' Result not transformed back to time domain')
end
