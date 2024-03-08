clc
clearvars

load('40cycles_avg4_20230222T093226_channel24.mat','RcvData','TX','doUseTxSweep','matchedFilterParams');
goodChannel = 1;
matchedFilterParams.activeRx = 24;
%%

bufferIndex = 1;

[totalSamples, numberOfChannels, numberOfFrames] = size(RcvData{bufferIndex});
numberOfTransmissions = length(TX);
numberOfTimeSamples = totalSamples ./ numberOfTransmissions;
rawData = reshape(RcvData{bufferIndex},[numberOfTimeSamples numberOfTransmissions numberOfChannels numberOfFrames]);
% rawData = squeeze(rawData);

% matched filtering
if doUseTxSweep
    bScan = applyMatchedFiltering(RcvData{bufferIndex}, matchedFilterParams);
end

% bScan dimensions: # time samples, channel, measurement number (at 50 fps)

%% axes
samplingFrequency = 31.25e6;
t=(0:size(rawData,1)-1)./samplingFrequency;
depth = t.*1500./2; % from two-way travel time in water to distance

channels = 1:32;
lateralPos = (channels-1)*2.7e-3; % 2.7 mm pitch

time = (0:size(bScan,3)-1) ./ 50; % time of measurement (with 50 fps)

%% plot evolvement of one channel over time to show motion

envelopeBestChannel = abs(hilbert(squeeze(bScan(:,goodChannel,:)))); % first measurement (all channels)
envelopeBestChannel = envelopeBestChannel ./ max(envelopeBestChannel(250:end,:),[],[1 2]);

%% Fig. 5d

figure(1)
imagesc(time,depth*1e3,20*log10(envelopeBestChannel)) % plot 
colormap bone
colorbar
ylim([0 30])
clim([-40 0])
xlabel('time [s]')
ylabel('depth [mm]')
title(['measurement over time, channel ' num2str(goodChannel)])

%% Fig. S13b

figure(2)
imagesc(time,depth*1e3,20*log10(envelopeBestChannel)) % plot 
colormap jet
colorbar
ylim([0 30])
clim([-40 0])
xlabel('time [s]')
ylabel('depth [mm]')
title(['measurement over time, channel ' num2str(goodChannel)])

%%

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
