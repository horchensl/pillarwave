
filenames = {'20cycles_avg8_ID2_0_40_60BPMall_20230324T115951.mat', ...
             '20cycles_avg8_ID2_00_45_60BPMall_20230324T120741.mat', ...
             '20cycles_avg8_ID2_20_60_60BPMall_20230324T120208.mat', ...
             '20cycles_avg8_ID2_30_100_60BPMall_20230324T121136.mat', ...
             '20cycles_avg8_ID2_75_130_60BPMall_20230324T121351.mat'};

labels = {'0 to 40 mmHg', ...
    '0 to 45 mmHg', ...
    '20 to 60 mmHg', ...
    '30 to 100 mmHg', ...
    '75 to 130 mmHg'
    };

goodChannels = [24 24 24 24 24];

Nfiles = length(filenames);

raw = []; % single channel time series, 10 seconds at 50 fps, four experiments

for fileIndex=1:Nfiles

    load(filenames{fileIndex},'RcvData','TX','doUseTxSweep','matchedFilterParams');
    
    %%
    
    bufferIndex = 1;
    
    [totalSamples, numberOfChannels, numberOfFrames] = size(RcvData{bufferIndex});
    numberOfTransmissions = length(TX);
    numberOfTimeSamples = totalSamples ./ numberOfTransmissions;
    rawData = reshape(RcvData{bufferIndex},[numberOfTimeSamples numberOfTransmissions numberOfChannels numberOfFrames]);
    rawData = squeeze(rawData);
    
    % matched filtering
    if doUseTxSweep
        bScan = applyMatchedFiltering(RcvData{bufferIndex}, matchedFilterParams);
    end
    
    % bScan dimensions: # time samples, channel, measurement number (at 50 fps)
    
    %% axes
    if fileIndex==1

        samplingFrequency = 31.25e6;
        t=(0:size(rawData,1)-1)./samplingFrequency;
        depth = t.*1540./2; % from two-way travel time in tissue to distance
        
        channels = 1:32;
        lateralPos = (channels-1)*2.7e-3; % 2.7 mm pitch
        
        time = (0:size(bScan,3)-1) ./ 50; % time of measurement (with 50 fps)
    end
    
    %% remove noise common to all channels
    
    bScan = bScan - mean(bScan(:,12:end,:),2);
    
   
    %% store 

    raw = cat(3,raw,squeeze(bScan(:,goodChannels(fileIndex),:)));

end

