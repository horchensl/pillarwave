clearvars
close all

for jobnr = 5 % loop over experiment number
    % default values, adjust when needed in the case-structure
    dx          = 75e-6;  % [m]
    dy          = 75e-6;  % [m]
    Cp          = 1484;   % [m/s]
    sens        = 47e-9;  % [nV/Pa]  sensitivity needle hydrophone = 47nV/Pa at 10MHz as default start value
    set_fMin    =  3.5e6; % lower frequency limit
    set_fMax    = 15.0e6; % upper frequency limit

    switch jobnr
        case 5 % array ID1, 1467x307 scan 2023_0117
            swpname     = 'Swp3-30MHz_Slope-5_T20.0us_Fs200MHz_universal.mat';
            downs_swp   = 2;                    % if 2 then downsample the sweep from 200MHz to 100MHz
            aveX        = 278 + (-5:5);         % area to average 'breedte'
            aveY        = 105 + (-5:5);         % 'hoogte'
            wstart      = 1700;                 % start windowing cor acou[sample nr] CAREFULL, Z-deviation in this measurement
            wlength     = 600;                  % window length of cor acou [number of samples]
            gain        = 10^(40/20);           % gain volt/volt of panametrics 5900R preamp
            n1st        = 307;                  % [number of samples] = figure Y-axis
            n2nd        = 1467;                 % [number of samples] = figure X-axis
            probe       = 10;                   % BNC probe attenuation
            do_sn3350   = false;                % sn2852 old hydrophone
            xlim_eff    = [3 30];
            max_freq    = 30e6;
            fnameacou   = 'fnameacous_1467X_307Y.mat';
            fnameexc	= 'fnameawg_1467X_307Y.mat';

    end

    %%  load AWG-sweep, Acoustics and Excitation data
    loaded      = load(swpname); % load chirp
    swp         = loaded.data;
    if downs_swp==2
        swp = swp(1:2:end);
    end
    load(fnameacou); % load raw data
    dat         = single(dat / (gain*sens)); % roughly correct for gain and sensitivity (later will be corrected for frequency dependence)
    dt			= t(2)-t(1); % define time step
    Fs			= 1/dt; % define sample frequency
    load(fnameexc,'exc'); % load excitation Voltage
    exc         = single(exc*probe); % Correct excitation Voltage for electrical probe setting

    %%  compress chirps on received and excitation waveforms
    swp         = swp(1:size(dat,1)); % force equal samples numbers in swo and acou-data
    cordat      = zeros(size(dat)); % init variable
    % correllate data with sweep
    for n=1:size(dat,2)
        cordat(:,n) = cor(dat(:,n),swp);
    end
    corexc      = zeros(size(exc)); %init variable
    % correlate recorded excitation signals with sweep
    for n=1:size(exc,2)
        corexc(:,n) = cor(exc(:,n),swp);
    end
    corexc = single(corexc); % make single precision
    cordat = single(cordat); % make single precision

    nt          = size(corexc,1); % number of time samples
    %dt          = t(2) - t(1);
    corexc      = ftshift(corexc,1,nt/2);                           % shift corexc in time to middle of scan-time
    tt          = t-(nt/2)*dt;                                      % create new time axis with 0 at middle

    %%  Windowing and FFT acoustics and excitation, reshape to cube
    L       = size(cordat,1);
    cordatw = wnd(cordat  ,4,wstart,wlength);
    corexcw = wnd(corexc  ,4,1759,600);

    dat     = reshape(dat,L,n1st,n2nd);
    exc     = reshape(exc,L,n1st,n2nd);
    cordat  = reshape(cordat,L,n1st,n2nd);
    corexc  = reshape(corexc,L,n1st,n2nd);
    cordatw = reshape(cordatw,L,n1st,n2nd);
    corexcw = reshape(corexcw,L,n1st,n2nd);

    %dt      = 1/Fs;
    f       = mkf(L,dt);                % make FFT freq vector
    f       = f(1:L/2);
    COR     = fft(cordatw);             % do FFT, capitals for f-domain
    EXC     = fft(corexcw);
    COR     = single(COR(1:L/2,:));     % throw-away negative freqs, make single to save memeory
    EXC     = single(EXC(1:L/2,:));
    COR     = reshape(COR,L/2,n1st,n2nd);
    EXC     = reshape(EXC,L/2,n1st,n2nd);

    %%  localize maximum in middle of selected area to pinpoint wavefield extrapolation distance dz
    [maxcoramp, midx]   = max(abs(cordatw(:,mean(aveY),mean(aveX))));
    dz                  = t(midx)*Cp;

    %% apply inverse wavefield-extrapolate , for all freqs
    wf_input = COR./EXC; % divide recorded data by excitation data
    PaV    = waveextr3d_fdom(wf_input,dx,dy,f,Cp,dz,'inverse',[0 0]); % !! normalized to EXC, both phase and magnitude

    %% correct for sensitivity curve of needle hydrophone and water absorbtion
    attenuation_dB      = 0.002.*((f'/1e6).^2).*dz*100; % calculate attenuation in fresh water @ 20 degrees C
    attenuation_lin     = 10.^(-attenuation_dB./20); % calculate linear attenuation vector
    correction          = sens./(Hydrophone2852_interp(f') .* attenuation_lin); % correct needle hydrophone sensitiivty curve for attenuation and frequency dependence sensitivity
    PaV                 =    PaV .*  correction'; % compensate recorded wavefield divided by excitation data, which was inverse wavefield extrapolated to location of aperture by needle hydrophone sensitivity and attenaution

    x                   = (0:size(PaV,3)-1) * dx;    % [m] PaV(f,y,x)
    y                   = (0:size(PaV,2)-1) * dy;    % [m]
    area2               = mean(abs(PaV(:,aveY,aveX)),2); %aveX and aveY are the pixel selection variables
    area3               = mean(area2,3);
    fmin                = find(f>set_fMin,1,'first'); % find index of lowest frequency
    fmax                = find(f>set_fMax,1,'first'); % find index of maximum freqeuncy
    [max_ave, maxidx]   = max(area3(fmin:fmax));        % maximum Pa/V averaged over area
    fmax_ave_idx        = maxidx + fmin -1;
    fmax_ave            = f(fmax_ave_idx);              % maximum occurs at this freq

    %% Make plots
    figure(1)
    subplot(311)
    plot(t*1e6,exc(:,mean(aveY),mean(aveX)))
    grid on; grid minor
    xlabel 'time [\mus]'
    ylabel 'Volt'
    title 'Excitation waveform (middle of selected area)'
    set(gca,'fontsize',20)
    subplot(312)
    plot(tt*1e6 ,corexc(:,mean(aveY),mean(aveX)))
    hold on
    plot(tt*1e6 ,corexcw(:,mean(aveY),mean(aveX)))
    grid on; grid minor
    xlabel 'time [\mus]'
    ylabel ''
    set(gca,'fontsize',20)
    subplot(313)
    plot(f/1e6,20*log10(abs(EXC(:,mean(aveY),mean(aveX)))))
    grid on; grid minor
    xlabel 'freq [MHz]'
    ylabel 'dB V'
    set(gca,'fontsize',20)

    figure(2)
    subplot(311)
    plot(t*1e6,dat(:,mean(aveY),mean(aveX)))
    grid on; grid minor
    xlabel 'time [\mus]'
    ylabel 'Pa'
    title 'Hydrophone waveform (middle of selected area)'
    set(gca,'fontsize',20)
    subplot(312)
    plot(t*1e6,cordat(:,mean(aveY),mean(aveX)))
    hold on
    plot(t*1e6,cordatw(:,mean(aveY),mean(aveX))) % trigger jitter corrected and windowed
    grid on; grid minor
    xlabel 'time [\mus]'
    ylim(maxcoramp* [-1.2 1.2])
    set(gca,'fontsize',20)
    subplot(313)
    plot(f/1e6,20*log10(abs(COR(:,mean(aveY),mean(aveX)))))
    grid on; grid minor
    xlabel 'freq [MHz]'
    ylabel 'dB Pa'
    set(gca,'fontsize',20)

    figure(3)
    imagesc(x*1e3,y*1e3,abs(squeeze(PaV(fmax_ave_idx,:,:))));
    title([num2str(dz*1e3,'%2.1f') ' mm'])
    colormap jet
    hC2 =colorbar;
    ylabel(hC2,'Pa/V')
    axis image
    clim('auto')
    set(gca,'fontsize',14)
end  % end jobnr-loop