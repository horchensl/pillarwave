clc
close all
clearvars

%%

load('Gaussian_pulse_20_cycles_MFapplied_20221109T123157_images.mat')

%%

interpolationFactor = 3;

centerIndex = round(length(xaxis)/2);

image = imageData(:,:,1);
maxAmplitude = max(max(image(15:end,centerIndex-15:centerIndex+15)));

fig = figure(1);
fig.Position = [100 100 600 600];

% perform cubic interpolation
[interpolatedImage,xInterpolated,zInterpolated] = interpolateImage(image,xaxis,zaxis,interpolationFactor);

imagesc(xInterpolated*1e3,zInterpolated*1e3,20.*log10((interpolatedImage)./maxAmplitude))
colormap gray
colorbar
axis equal tight
clim([-40 0])
xlabel('lateral position [mm]')
ylabel('depth [mm]')
ylim([0 50])
title('Phantom image')

