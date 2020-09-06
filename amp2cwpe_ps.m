function [cwpe, fa] = amp2cwpe_ps(B1max,N_pulses,pulse_duration,pulse_interval,spoiler)
%[cwpe, fa] = amp2cwpe_ps(B1max,N_pulses,pulse_duration,pulse_interval,spoiler)
%
% Input:
% B1max - B1 amplitude [uT]
% N_pulses - number of pulses in the pulse train
% pulse_duration - single pulse duration [ms]
% pulse_interval - interpulse delay [ms]
% spoiler - final spoiler (before readout) [ms]
%
% Output:
% cwpe - B1 average [uT]
% fa - flip angle
%
% Example:
% [cwpe, fa] = amp2cwpe_ps(2,20,50/1000,50/1000,4/1000)
%

% 20191106
% (c) Vitaliy Khlebnikov, PhD (Radboud University)

pw=pulse_duration;
tr=pulse_duration+pulse_interval;
tr2=pulse_duration+pulse_interval+spoiler;

saturation = [0, 413, 859, 1337, 1847, 2389, 2962, 3568, 4204, 4870,...
    5565, 6288, 7038, 7814, 8613, 9435, 10277, 11138, 12015, 12906,...
    13809, 14721, 15640, 16563, 17487, 18409, 19327, 20237, 21136, 22023,...
    22892, 23742, 24570, 25372, 26146, 26889, 27598, 28271, 28905, 29497,...
    30046, 30549, 31005, 31411, 31767, 32070, 32320, 32515, 32655, 32739,...
    32767, 32739, 32655, 32515, 32320, 32070, 31767, 31411, 31005, 30549,...
    30046, 29497, 28905, 28271, 27598, 26889, 26146, 25372, 24570, 23742,...
    22892, 22023, 21136, 20237, 19327, 18409, 17487, 16563, 15640, 14721,...
    13809, 12906, 12015, 11138, 10277, 9435, 8613, 7814, 7038, 6288,...
    5565, 4870, 4204, 3568, 2962, 2389, 1847, 1337, 859, 413, 0 ]; %Sinc-pulse

%saturation=ones(size(saturation)); % for Rectangular pulses

    rf=saturation;
    rf = rf/max(rf);
    dt = pw/length(rf);

    B1max=B1max/100;    
    cwpe_1=B1max.*(sqrt(sum(rf.^2).*dt.*4257.7^2/tr)/42.577);
    fa=B1max.*(4257*360*sum(rf).*dt);  %in Guass
    
    cwpe_2=B1max.*(sqrt(sum(rf.^2).*dt.*4257.7^2/tr2)/42.577);
    
    cwpe=((N_pulses-1)*cwpe_1+cwpe_2)/N_pulses;

end