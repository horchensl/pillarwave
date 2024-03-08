function output = applyMatchedFiltering(RcvBuffer, matchedFilterParams)
    nt = matchedFilterParams.nt;
    nf = matchedFilterParams.nf;
    S0 = matchedFilterParams.S0;
    activeRx = matchedFilterParams.activeRx;
    nSamples = matchedFilterParams.nSamples;
    nTransmissions = matchedFilterParams.nTransmissions;

    numberOfFrames = size(RcvBuffer,3);

    nActiveChannels = length(activeRx);
    % put the channels from all transmissions side-by-side
    dataRx = reshape(RcvBuffer(:,:),[nt nTransmissions*nActiveChannels*numberOfFrames]);
	dataRx(1:10,:) = 0; % set intial pulse to zero
    Rf = fft(dataRx,nf);
    Rf = bsxfun(@times,Rf,S0); % apply the matched filter
    compressed = real(ifft(Rf));
    compressed = compressed(1:nt,:); % cut off excess samples

    output = reshape((compressed),[nSamples,nActiveChannels,numberOfFrames]); % place back into original array

end