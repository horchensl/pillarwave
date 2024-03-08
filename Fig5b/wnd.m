function y = wnd(x,w_type,w_start,w_length)
%WND	WND(x,type,start,length) windows or tapers signal x over
%	'length' points, starting at point 'start'.
%TSMAT
%TLBX	If 'type'is positive the windowing is done in the time
%	domain, if 'type' is negative in the frequency domain.
%       In the latter case, positive and negative frequencies are
%	treated equally, thus preserving the realness of the time
%	domain signal. Also: 'start' + 'length' <= N/2 + 1 !!!
%
%	The signal is returned in the original domain.
%
%	The following window types are possible:
%
%	type = 1, symmetric cosine window over 'length', data before
%		  'start' and after 'start + length' is zeroed.
%	type = 2, left cosine window over 'length', data before 'start'
%		  is zeroed, data after 'start + length' is not affected.
%	type = 3, right cosine window over 'length', data before 'start'
%		  is not affected, data after 'start + length' is zeroed.
%	type = 4, flat top window with cosine flanks of 'length/10' length,
%		  data before 'start' and after 'end' is zeroed.
%   type = 5, same as 4 but cosine flanks of 'length/20'.
%
%	NB. Applied columnwise only !!!
%
%	See also FILTER.

%	Uilke Stelwagen, September 1991, June 1996.
%	Copyright (C) 1991, 1996, Institute of Applied Physics, TNO-TPD,
%	The Netherlands.

if (nargin < 4)
   help wnd
   return
end

% Check & prepare

fw_bw = 0;
[nr nc] = size(x);
l_x   = nr;
w_end = w_start + w_length - 1;

if(w_start<1)
  disp('  Startpoint redefined to 1 !')
  disp('  Windowlength appropriately redefined !')
  w_start = 1;
  w_length = w_end - w_start + 1;	% Don't forget windowlength !
end

if(w_type>0)
  if(w_end>l_x)
    disp('  Endpoint redefined to lenght(input) !')
    disp('  Windowlength appropriately redefined !')
    w_end = l_x;
    w_length = w_end - w_start + 1;	% Don't forget windowlength !
  end
elseif(w_type<0)		% Frequency windowing, take care !
  if(w_end>(l_x/2))
    disp('  Frequency windowing !')
    disp('  Endpoint redefined to lenght(x)/2 !')
    w_end    = l_x/2;
    w_length = w_end - w_start + 1;	% Don't forget windowlength !
  end
end

if( (w_type>0) & (iscmplx(x)) )
  x = real(ifft(x));
  fw_bw = 1;			% Remember to forward FFT the result,
elseif( (w_type<0) & (isreal(x)) )
  x = fft(x);
  fw_bw = -1;			% or to backward FFT the result !
end

% Build window

w = ones(nr,1);	% Note that here w is 1 column wide !!!

if(w_type==1 | w_type==-1)	% Symmetric cosine window
  if(w_start>1)
    w(1:w_start-1) = zeros(w_start-1,1);
  end
  w(w_start:w_end) = .5  - .5*cos(2*pi*(0:w_length-1)/(w_length-1));
  if(w_end<l_x)
    w(w_end+1:l_x) = zeros(l_x-w_end,1);
  end
elseif(w_type==2 | w_type==-2)	% Left cosine window
  if(w_start>1)
    w(1:w_start-1) = zeros(w_start-1,1);
  end
  w(w_start:w_end) = .5  - .5*cos(pi*(0:w_length-1)/(w_length-1));
elseif(w_type==3 | w_type==-3)	% Right cosine window
  if(w_end<l_x)
    w(w_end+1:l_x) = zeros(l_x-w_end,1);
  end
  w(w_start:w_end) = .5  + .5*cos(pi*(0:w_length-1)/(w_length-1));
elseif(w_type==4 | w_type==-4)	% Flat top cosine flanked 1,8,1
  f_length = fix(w_length/10);
  if(w_start>1)
    w(1:w_start-1) = zeros(w_start-1,1);
  end
  w(w_start:(w_start+f_length-1)) = .5  - .5*cos(pi*(0:f_length-1)/(f_length-1));
  w((w_end-f_length+1):w_end) = .5  + .5*cos(pi*(0:f_length-1)/(f_length-1));
  if(w_end<l_x)
    w(w_end+1:l_x) = zeros(l_x-w_end,1);
  end
elseif(w_type==5 | w_type==-5)	% Flat top cosine flanked 0.5,9,0.5
  f_length = fix(w_length/20);
  if(w_start>1)
    w(1:w_start-1) = zeros(w_start-1,1);
  end
  eps=1e-6;
  w(w_start:(w_start+f_length-1)) = .5  - .5*cos(pi*(0:f_length-1)/(f_length-1+eps));
  w((w_end-f_length+1):w_end) = .5  + .5*cos(pi*(0:f_length-1)/(f_length-1+eps));
  if(w_end<l_x)
    w(w_end+1:l_x) = zeros(l_x-w_end,1);
  end
else
  disp('  Sorry, unknown window type. Try help WND.')
  return
end

% If necessary adjust for frequency windowing

if(w_type<0)
  w((l_x/2+1):1:l_x) = w(l_x/2:-1:1); % This mirrors the positive
end                                   % freq. onto the negative !

% Duplicate columns

w = w(:,ones(1,nc));	% Tony's trick !

% Apply window and if necessary transform

y = x.*w;

if(fw_bw==1)
  y = fft(y);
elseif(fw_bw==-1)
  y = real(ifft(y));
end
