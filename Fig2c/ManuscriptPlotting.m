% Matlab script to generate Manuscript figure 2c using RAW data
%
% Files needed:
% - 'TXwaveform2c.mat'     :  Transmit wavefrom function [-]
% - 'TxData2c.mat'         :  Measured transmit waves at transducer terminals,
%                             compensated for 1:10 high impedance probe
%                             [Volt], including time traces
% - 'HydrophoneData2c.mat' :  Measured hydophone response [Volt], including time traces
%
% 7 March 2024
% TNO, Oude Waalsdorperweg 63, Netherlands
% Copyright (C) 2024 Institute of Applied Physics, TNO, The Netherlands.

clearvars

%% Known acquisition parameters and processing settings
fsample = 100e6;                    % Acquisition sample frequency
Receive_amplifier_gain_dB = 40;     % Received signal gain in [dB]
Hydrophone_sens = 4.5937e-08;       % Hydrophone sensitivity @ 8.2 MHz [V/Pa]
Hydrophone_Distance = 0.0298;       % Hydophone distance [m] from transdcuer surface, derived from time of flight.
Freq_sel = 8.93e6;                  % Selected frequency [Hz]
c0 = 1480;                          % Speed of sound [m/s]

dx = 75e-6;                         % Scan step  X-axis [m]
dy = 75e-6;                         % Scan step  Y-axis [m]
Nsteps_X = 281;                     % Number of scan positions, X-axis
Nsteps_Y = 241;                     % Number of scan positions, Y-axis
x = [0:Nsteps_X - 1].*dx - 10.5e-3; % Absolute scan positions X-axis [m]
y = [0:Nsteps_Y - 1].*dy - 9e-3;    % Absolute scan positions y-axis [m]
X_loc_offset = 1.25;                % Centre alignment offset [mm], X-axis
Y_loc_offset = 0;                   % Centre alignment offset [mm], Y-axis

NFFT = 4096;            % Number of points used for fourier transform
N_window = NFFT/8;      % Number of points used for zero padding and windowing time trace
range0 = 1984:3007;     % Range of sample points containing pressure wave (i.e. time windowing time trace)  

%% Load Raw data
load('TXwaveform2c.mat');               % Lineair frequency sweep transmit waveform 
load('TxData2c.mat');                   % Transmit signal on transducer connector, compensated for 1:10 probe attenuation factor [Volt] and time traces [sec]
us_src_sweep = mean(d_allsweep_src,2);  % Averaging transmit signal across all scan positions, assuming stable excitation

load('HydrophoneData2c.mat');           % Hydophone output [Volt] and time traces [sec]
d_allsweep_us = 1/Hydrophone_sens.*10^(-Receive_amplifier_gain_dB./20).*d_allsweep_us; % Pressure conversion [Volt]-> [Pa]

%% Collapse by correlating with transmit waveform 
cor_resp = zeros(N_window+min(size(d_allsweep_us)),max(size(d_allsweep_us))); % Memory allocation

for i = 1:max(size(d_allsweep_us))
    temp = xcorr([zeros(1,N_window) d_allsweep_us(:,i).'].',SweepData);
    cor_resp(:,i) = temp(length(d_allsweep_us(:,i).')+N_window:end);
end

cor_src_sweep = 0.*cor_resp; % Memory allocation
for i = 1:max(size(d_allsweep_src))
    temp = xcorr([zeros(1,N_window) d_allsweep_src(:,i).'].',SweepData);
    cor_src_sweep(:,i) = temp(length(d_allsweep_src(:,i).')+N_window:end);
end

%% Transfer function and compensations
index = round(Freq_sel*NFFT/fsample) + 1; % Selecting the sample representing 8.2 MHz

Corrected =  reshape(cor_resp, 4608, Nsteps_Y , Nsteps_X);

FFT_resp = fft(Corrected(range0,:,:),NFFT);         % Spectrum of hydrophone response data
FFT_src = fft(cor_src_sweep(1:NFFT/4,:),NFFT);      % Spectrum of source data
FFT_src = reshape(FFT_src,NFFT,Nsteps_Y,Nsteps_X);
FFT_transfer = squeeze(FFT_resp(index,:,:))./squeeze(FFT_src(index,:,:)); % Calculate transfer function [Pa/Volt]

% Frequency dependent attenuation in water :: = 0.002 dB/MHz^2/cm
attenuation_dB = 0.002.*((Freq_sel/1e6).^2).*Hydrophone_Distance*100;
attenuation_lin = 10.^(-attenuation_dB./20); 
FFT_transfer = FFT_transfer./attenuation_lin .'; % Correct for attenuation [Pa/Volt]

%% Wavefield Extrapolation
N = max(Nsteps_X, Nsteps_Y); % Maximum grid size
N = 2^(round(log2(N))+1);

range_X_dat2 = N/2-floor(Nsteps_X/2):N/2-floor(Nsteps_X/2)+Nsteps_X-1;  % Make grid
range_Y_dat2 = N/2-floor(Nsteps_Y/2):N/2-floor(Nsteps_Y/2)+Nsteps_Y-1;  % Make grid

dat_noise = abs(FFT_transfer([1 Nsteps_Y],1:Nsteps_X)); % Get noise data
dat = FFT_transfer(1:Nsteps_Y,1:Nsteps_X);              % Get transfer data
dat2 = mean(dat_noise(:)).*(randn(N) + randn(N).*sqrt(-1))*0.5*sqrt(2); % Make complex noise data with max level equal to mean value
dat2(range_Y_dat2 ,range_X_dat2) = dat;                 % Fill remaining area with random noise to remove edge effects

dkx = 1/(dx * N);
kx  = 2.*pi.*[(0:N/2-1) ((N/2:N-1)-N)]'*dkx;    % Spatial frequency range
[kx,ky]=meshgrid(kx,kx);                        % Make mesh grid

w=2*pi*Freq_sel;     % radial frequency
dat2=fft2(dat2);     % 2D spatial Fourier transform

W = sqrt( (w./c0).^2 - kx.^2 - ky.^2);              %
W = sqrt(-1).*abs(real(W)).*sign(w) - abs(imag(W)); % Rayleigh II operator in kx-w domain
dat2=exp(W.*Hydrophone_Distance).*dat2;             % Inverse wave extrapolation
dat_extrap=ifft2(dat2);                             % Inverse 2D spatial Fourier transform


%% Plotting 
imagesc(([0:N-1]-N/2).*dx*1000-X_loc_offset,([0:N-1]-N/2).*dy*1000-Y_loc_offset,abs(dat_extrap)), hold on
axis( 'equal' );
axis xy
clim([0 1500])
xlabel( 'Lateral axis [mm]','FontSize',11 );
ylabel( 'Elevation axis [mm]','FontSize',11 );
xlim([-1 1].*10.5)
ylim([-1 1].*10.5)
c = colorbar;
c.Label.String =  ['Transmit Sensitivity ' '[Pa/V]'];
colormap jet,
label = { c.Label.String,sprintf('Extrapolated field @ %.2f MHz' ,Freq_sel/1e6)};
title(label,'FontSize',11)
drawnow