function HydSens_int = Hydrophone2852_interp(Freq_target)
% Gets hydrophone sensitivity at custom frequencies. Uses interp1 to
% interpolate the calbration data.

load('Hydrophone_sn2852.mat');

% moving average filter
temp = Hydrophone_sn2852.Sensitivity.SensV_Pa;
N = 3;
for j=1:2,
    for i = 1:length(temp),
        if i<=N/2,
            temp(i) = mean(temp(i:i+floor(N/2)));
        elseif i >= length(temp)-N/2,
            temp(i) = mean(temp(i-floor(N/2):i));
        else
            temp(i) = mean(temp(i-floor(N/2):i+floor(N/2)));
        end
    end
    temp = flipud(temp);
end

% Extend curve to fit target frequency range
df = 1e6*mean(diff(Hydrophone_sn2852.Sensitivity.FreqMHz));
Diff = diff(temp);
rc_start = mean(Diff(1:5))./df;
rc_end = mean(Diff(end-5:end))./df;


F_0 = Hydrophone_sn2852.Sensitivity.FreqMHz(1).*1e6;
F_end = Hydrophone_sn2852.Sensitivity.FreqMHz(end).*1e6;
if min(Freq_target) < F_0,  %#ok<*BDSCI>
    f0 = min(Freq_target);
    Ftemp0 = fliplr(F_0:-df:f0-df);
    temp0 = rc_start.*(Ftemp0(1:end-1)-F_0)+temp(1); 
    temp = [temp0 temp.'];
    Ftemp =  [Ftemp0(1:end-1) Hydrophone_sn2852.Sensitivity.FreqMHz.'*1e6];
else
    Ftemp =  Hydrophone_sn2852.Sensitivity.FreqMHz.'*1e6;
end

if max(Freq_target)>F_end, 
    f_end= max(Freq_target);
    Ftemp0 = F_end:df:f_end+df;
    temp0 = rc_end.*(Ftemp0(2:end)-F_end)+temp(end);
    [n,m]=size(temp);
    if n>m,
        temp = temp.';
    end
    temp = [temp temp0];
    Ftemp =  [Ftemp Ftemp0(2:end)];
end

% Interpolate curve to specific frequencies
HydSens_int = interp1(Ftemp, temp, Freq_target);



    