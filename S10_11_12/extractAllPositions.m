clc

window = zeros(4,4);
window(:,1) = [12.5 12.9 13.2 13.8]*1e-3; % depth
window(:,2) = [14.5 14.8 15.5 15.6]*1e-3; % depth
window(:,3) = [20.5 21.0 21.8 22.3]*1e-3; % depth
window(:,4) = [23.0 23.4 24.2 24.6]*1e-3; % depth
Nsamples = 500;

wallDepths = zeros(Nsamples,4,Nfiles);


for depthIndex = 1:4
    windowFunction = squaredcos(depth,window(:,depthIndex)).';
    referenceTrace=raw(:,1,1).*windowFunction; % relative to first trace of first file
    for fileIndex=1:Nfiles
        
        part = raw(:,1:Nsamples,fileIndex).*windowFunction;
        
        [shift,envInterpolated,maxIndex,A] = crossCorrelation(part,referenceTrace);
        
        ddepth = depth(2)-depth(1);
        
        wallDepths(:,depthIndex,fileIndex) = ddepth.*shift;
        
    end

end


%% selection of 5 seconds

Nsamples = 250;
startIndices = [5 158 117 61 193];

innerDiameter = zeros(Nsamples,Nfiles);
outerDiameter = zeros(Nsamples,Nfiles);
for n=1:Nfiles
    innerDiameter(:,n) = wallDepths(startIndices(n):startIndices(n)+Nsamples-1,3,n) - wallDepths(startIndices(n):startIndices(n)+Nsamples-1,2,n);
    outerDiameter(:,n) = wallDepths(startIndices(n):startIndices(n)+Nsamples-1,4,n) - wallDepths(startIndices(n):startIndices(n)+Nsamples-1,1,n);
end


