%wave extrapolation using Rayleigh II operator in kx-w domain.
%
% input: dat    : data in f,x,y domain (even number of spatial samples)
%        dx     : spatial sampling interval
%        dy     : spatial sampling interval
%        f      : frequencies
%        c      : medium velocity
%        dz     : extrapolation distance
%        type   : extrapolation direction inverse/forward
%        Ntap   : taper points in [x y ang0 dang cmin cmax] (optional)
%
%
% usage
%  [dat]=waveextr3d_fdom(dat,dx,dy,f,c,dz,type,Ntap);
%
% type=ínverse/forward
function [dat]=waveextr3d_fdom(dat,dx,dy,f,c,dz,type,Ntap);



% fs=1/dt;
[N,M,O]=size(dat);

% if rem(N,2)~=0,
%   dat=dat(1:N-1,:,:);
%   N=N-1;
% end
if rem(M,2)~=0,
    dat=dat(:,1:M-1,:);
    M=M-1;
end
if rem(O,2)~=0,
    dat=dat(:,:,1:O-1);
    O=O-1;
end

if nargin<8
    win=1;
    Ntap=[];
else
    winx=mktaper(M,Ntap(1));
    winy=mktaper(O,Ntap(2));
    win=winx(:)*winy(:).';
    for l=1:N
        dat(l,:,:)=squeeze(dat(l,:,:)).*win;
    end
end


X=(0:M-1).*dx;
Y=(0:O-1).*dy;
% T=(0:N-1).*dt;

% apply angle/velocity filter
if length(Ntap)==6
    [win]=kxky_angle_filter(mkx(O,dy),mkx(M,dx),Ntap(3),Ntap(4),f,[Ntap(5) Ntap(6)],0);
else
    win=1;
end
kx=2.*pi.*mkf(M,dx);
kx=kx*ones(1,O);
kx=kx;

ky=2.*pi.*mkf(O,dy);
ky=ky*ones(1,M);
ky=ky';

w=2*pi.*f;
% w=w*ones(1,M);

% forward FT
dat=fft(dat,[],2);
dat=fft(dat,[],3);

dat=dat.*win;

for l=1:length(w),
    W=sqrt( (w(l)./c).^2 - kx.^2 - ky.^2) ;
    W=j.*abs(real(W)).*sign(w(l)) - abs(imag(W));
    
    if strcmp(type,'inverse'),
        dat(l,:,:)=exp(W.*dz).*squeeze(dat(l,:,:));
    else
        dat(l,:,:)=exp(conj(W).*dz).*squeeze(dat(l,:,:));
    end;
    
end
% inverse FT
dat=ifft(dat,[],3);
dat=ifft(dat,[],2);
% dat=real(ifftn(dat));
