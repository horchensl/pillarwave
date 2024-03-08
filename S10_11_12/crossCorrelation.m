function [shift,envInterpolated,maxIndex,A] = crossCorrelation(A,referenceTrace, filterFreqs, filterOrder)

if nargin < 4
    filterOrder = 4;
end

if nargin > 2
    [b,a]   = butter(filterOrder,filterFreqs);
    A     = filtfilt(b,a,A);
end
%%


xc = zeros(2*size(A,1)-1,size(A,2));

for step=1:size(A,2)
    cross = xcorr(A(:,step),referenceTrace,'none'); % use fixed trace as reference
    xc(:,step) = cross;
end


%% determine envelope and normalize per frame 

env = abs(hilbert(xc));
env = xc./max(env,[],1);

%% interpolate to find sub-sample peak

Nc = floor(size(env,1)/2); % number of samples in the correlation
t=-Nc:Nc;

interpolationFactor=10;
[tInterpolated,envInterpolated] = sincInterpolationFourier(t,env,interpolationFactor);

[~,maxIndex] = max(envInterpolated);
shift = tInterpolated(maxIndex);

%% reference peak

[~,maxPositionInReference] = max(abs(hilbert(referenceTrace)));
shift = shift+maxPositionInReference;

