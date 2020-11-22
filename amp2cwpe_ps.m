function [cwpe, fa] = amp2cwpe_ps(mt_rf_pulse_profile, b1_amp_mt, N_mt_pulses, mt_pulse_dur, mt_pulse_int, final_spoiler)
%[cwpe, fa] = amp2cwpe_ps(mt_rf_pulse_profile, b1_amp_mt, N_mt_pulses, mt_pulse_dur, mt_pulse_int, final_spoiler)
%
% Script for pseudosteady state (ps) CEST-MRI sequences for conversion between b1+ field of MT
% pulses [amplitude in uT] and cwpe approximation
%
% Input:
% mt_rf_pulse_profile - mt rf pulse profile
% b1_amp_mt - b1+ field of MT pulses [uT]
% N_mt_pulses - number of MT pulses
% mt_pulse_dur - duration of MT pulses [s]
% mt_pulse_int - interval of MT pulses [s]
% final_spoiler - final spoiler duration [s] after MT preparation
%
% Output:
% cwpe - b1+ field cwpe
% fa - flip angle

pw=mt_pulse_dur;
tr=mt_pulse_dur+mt_pulse_int;
tr2=mt_pulse_dur+mt_pulse_int+final_spoiler;

rf=mt_rf_pulse_profile;
rf = rf/max(rf);
dt = pw/length(rf);

b1_amp_mt=b1_amp_mt/100;
cwpe_1=b1_amp_mt.*(sqrt(sum(rf.^2).*dt.*4257.7^2/tr)/42.577);
cwpe_2=b1_amp_mt.*(sqrt(sum(rf.^2).*dt.*4257.7^2/tr2)/42.577);

fa=b1_amp_mt.*(4257*360*sum(rf).*dt);  
cwpe=((N_mt_pulses-1)*cwpe_1+cwpe_2)/N_mt_pulses;

end