
% make window in kx-ky domain to select angle-wavenumber range
%
% input: kx: kx-axis created using mkkx.m
%        ky: ky-axis created using mkkx.m
%        ang0 : central angle to pass
%        dang : angle range to pass
%        f    : frequency
%        crng : cut-off velocity range
%
% usage: [win,winANG,winR]=kxky_angle_filter(kx,ky,ang0,dang,f,cmin,verbose);
%
%--------------------------------------------------------------------------

% by Arno Volker
% 28-01-2022

function [win,winANG,winR]=kxky_angle_filter(kx,ky,ang0,dang,f,crng,tukeyPar,verbose);


% make 2D grid
[KX,KY]=meshgrid((kx),(ky));

% angle and radius
ANG=atan2d(KY,KX);
KR=sqrt(KX.^2+KY.^2);

% angular window
ang_i=-180:.1:180;
Amp_i=cosd(ang_i)>=cosd(dang);
m=find(Amp_i==1);
win_i=zeros(size(Amp_i));
if length(m)<length(win_i)
    win_i(m)=tukeywin2(length(m),tukeyPar);
else
    win_i(m)=1;
end
AA=ANG-ang0;
AA=acosd(cosd(AA));
winANG=interp1(ang_i,win_i,AA);

% radial window
kr=linspace(-max(KR(:)),max(KR(:)),1e3);

switch(length(crng))
    case(1)
        cmin=crng;
        m=find(abs(kr)<=min([max(abs(kx)),max(abs(ky)),f./cmin]));
        1
    case(2)
        cmin=min(crng);
        cmax=max(crng);
        m=find(kr>=max(0,f./cmax) & (kr)<=min([max(abs(kx)),max(abs(ky)),f./cmin]));
end

win_r=zeros(size(kr));
win_r(m)=tukeywin2(length(m),tukeyPar);
winR=interp1(kr,win_r,KR);

win=winANG.*winR;

if verbose>0
    figure(11); clf
    subplot(221)
    imagesc(kx,ky,cosd(ANG-ang0)>=cosd(dang));  axis image
    xlabel('kx [m^{-1}]');  xlabel('ky [m^{-1}]');
    subplot(222)
    plot(ang_i-ang0,Amp_i)
    xlabel('angle [deg]');
    subplot(223)
    plot(ang_i-ang0,win_i)
    xlabel('angle [deg]');
    subplot(224)
    imagesc(fftshift(kx),fftshift(ky),fftshift(win)); axis image
    xlabel('kx [m^{-1}]');  xlabel('ky [m^{-1}]');
    colormap(jet)
    drawnow;
    
end

