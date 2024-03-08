function [tInterpolated,interpolatedSignal] = sincInterpolationFourier(t,signal,interpolationFactor)
S = fft(signal);
N=size(signal,1);
positiveN = floor(N/2)+1; % size of positive half of the spectrum

% pad with zeros in Fourier domain
if size(S,2) == 1
    S=interpolationFactor*[S(1:positiveN); zeros((interpolationFactor-1).*N,1); S(positiveN+1:end)];
else
    S=interpolationFactor*[S(1:positiveN,:); zeros((interpolationFactor-1).*N,size(S,2)); S(positiveN+1:end,:)];
end

interpolatedSignal = real(ifft(S));
dt = t(2)-t(1);
tInterpolated=(0:(interpolationFactor*N-1)).*dt/interpolationFactor + t(1); % time axis after interpolation
end