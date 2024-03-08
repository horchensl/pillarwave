function [interpolatedImage,xInterpolated,zInterpolated] = interpolateImage(image,xaxis,zaxis,factor)

% use cubic interpolation for upsampling an image by a given factor,
% return the upsampled axes as well

dx = xaxis(2)-xaxis(1);
dz = zaxis(2)-zaxis(1);

xInterpolated = xaxis(1):dx/factor:xaxis(end);
zInterpolated = zaxis(1):dz/factor:zaxis(end);

[X,Z] = meshgrid(xaxis,zaxis);
[Xout,Zout] = meshgrid(xInterpolated,zInterpolated);

interpolatedImage = interp2(X,Z,image,Xout,Zout,'cubic');
interpolatedImage = max(interpolatedImage,0); % assure positive image

end