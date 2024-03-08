%Apply time shift in the frequency domain
% input dat: data in time domain
%       dt : sampling interval
%       tau: time shift
%
%usage
%   dat=ftshift(dat,dt,tau);
%
%--------------------------------------------------------------------------
function  dat=ftshift(dat,dt,t);

inp_cmplx=iscmplx(dat); %check if input is complex

[m,n]=size(dat);


f=mkf(m,dt);
%w=2*pi*f*ones(1,n);
w=2*pi*f;


[m1,n1]=size(t);
if m1>1 & n1==1
   t=t.';
end
if  m1==1 & n1==1
    t=t*ones(1,n);
end

if inp_cmplx==false % real input must give real output
    dat=real(ifft(fft(dat).*exp(-j.*w*t)));
else
    dat=(ifft(fft(dat).*exp(-j.*w*t)));
end
