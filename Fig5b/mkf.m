function [f,df] = mkf(ns,dt)
%MKF	MKF(ns,dt) makes a vector f for use as frequency axis with TSMAT
%	plot routines. Optionally the frequency sample interval df is
%TSMAT	also given.
%TLBX	
%	Usage:	f = MKF(ns,dt);		or:	[f,df] = MKF(ns,dt);	
%
%	The number of samples is denoted by ns, while dt and df are
%	the temporal and frequency sample intervals respectively.
%
%	If ns is even, the axis runs from DC (f = 0, sample 1) up to
%  F_Nyquist-df (sample ns/2), and next from -F_Nyquist (sample ns/2 + 1)
%	to -df (sample ns), where F_Nyquist = 1/(2*dt).
%  For ns is odd, the axis goes from DC to F_Nyquist (sample
%  round(ns/2+eps) and next from -F_Nyquist (sample round(ns/2+eps) + 1)
%  to -df (sample ns).
%
%	PLS(f,r) or PLS(f,R) then plots the frequency spectrum of
%	vector r or R.
%
%	See also MKT, TSFMT, TSCRT and PLS.

%	Uilke Stelwagen, July 1991, revised June 1996 to comply with FFTSHIFT
%  and again in August 1999 to account for odd ns.
%	Copyright (C) 1991,1996,1999 Institute of Applied Physics, TNO-TPD,
%	The Netherlands.

if (nargin ~= 2)
   help mkf
   return
end

% Make a frequency axis with bins i according to:
%
%	f = (i - 1)*Fs/ns,		for i = 1 up to i = ns/2,
% and
%	f = (i - ns - 1)*Fs/ns,		for i = ns/2 + 1 to i = ns.
%
% Note special treatment of negative frequencies !

Fs = 1/dt;
df = Fs/ns;
%f  = [(0:ns/2-1) ((ns/2:ns-1)-ns)]'*df;	
f  = [(0:round(ns/2+eps)-1) ((round(ns/2+eps):ns-1)-ns)]'*df;	
