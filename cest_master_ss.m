function cest_master_ss
% Master script for modeling of CEST-MRI effects for steady state (ss)
% CEST-MRI sequences

% See Khlebnikov, V. et al. Comparison of pulsed three-dimensional
% CEST acquisition schemes at 7 tesla: steady state versus pseudosteady state.
% Magn. Reson. Med. 77, 2280–2287 (2017).

%
% (c) 2020, Vitaliy Khlebnikov, PhD
% 

close all
b_mz_evolution = 'false'; % SLOW! Show the signal evolution of Mz component. Options: false (for NO) or true (for YES)

% pools = {'Pool Name', Water T1 [s], Water T2 [s], Chemical Shift of the Pool relative to Water [ppm], ...
% Concentration of Exchangeable Protons [mM], Exchange Rate [Hz], 'LineShape'};

% Options for LineShape: 'lorentzian', 'gaussian' or 'super_lorentzian'

% The system is easily scalable, any number of pools can be added using
% the template above (each pool is a cell element)
pools={...
    {'H2O', 2.0,  40e-3,    0, 2*55.6*1000,  0, 'lorentzian'}...
    {'APT', 1.0,  10e-3,  3.5,       166.8, 50, 'lorentzian'}...
    {'MT',  1.0,  10e-6,   -5,       12232, 50, 'lorentzian'}...
    {'NOE', 1.0, 0.3e-3, -3,5,        6672, 10, 'lorentzian'}...
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CEST-MRI pulse sequence parameters
cest_sequence_type = 'ss'; % with EPI readout
main_field_strength = 7.0; % the main magnetic field strength [T]
% Sampling strategy explained:
% jump: 0, -1, 1, -2, 2, -5, 5 ...
% linear: -5, -2, -1, 0, 1, 2, 5 ...
% random: 2, 0, -5, 1, 5, -2 ...
sampling_strategy = 'linear'; % frequency sampling strategy. Options: jump, linear or random.
freq_offsets=[-5:1:5]; % offsets in ppm, slave script will generate the desired offsets depending on the sampling strategy specified
b1_amp_mt = 2; % b1+ field of MT pulses [uT]
mt_pulse_dur = 25/1000; % duration of MT pulses [s]
mt_pulse_int = 0; % interval of MT pulses [s]
final_spoiler = 0.5/1000; % final spoiler duration [s] after MT preparation
seq_tr = 65/1000; % sequence TR [s]
k0_time = 8.8; % the moment kspace center is reached [s]
ds_time = 16.9; % time per single offset acquisition [s]
N_seq_repeatitions = ceil(ds_time/seq_tr); % number of sequence repetitions (1 MT pulse per each sequence repetition)
N_dummy_offsets = 0; % number of dummy offsets acquisitions
t1_recovery = 10; %  post offset T1 recovery [s]
spacing_for_mz = ds_time + t1_recovery; % spacing out Mz progress (used only if b_mz_evolution = 'true')
mt_rf_pulse_file = 'am_sg_100_100_0.txt'; % % pulse shape of mt RF pulse
file_id = fopen(mt_rf_pulse_file,'r');
mt_rf_pulse_profile = fscanf(file_id,'%d'); % mt rf pulse profile
fclose(file_id);
%% Readout parameters (important to take into account if t1_recovery is short)
b1_amp_ex = 1.5; % b1+ field of readout excitation pulses [uT]
N_ex = 1; % number of readout excitations per shot
N_shots = 1; % number of shots
ex_dur = 0.2240/1000; % excitation duration [s]
read_dur = 8.3/1000; % EPI readout duration [s]
ex_rf_pulse_file = 'am_sg_400_150_125.txt'; % pulse shape of excitation RF pulse
file_id = fopen(ex_rf_pulse_file,'r');
ex_rf_pulse_profile = fscanf(file_id,'%d'); % excitation rf pulse profile
fclose(file_id);
%% Additional parameters
other_seq_par_container{1} = k0_time;
other_seq_par_container{2} = ds_time;
other_seq_par_container{3} = mt_rf_pulse_profile;
seq_dead_time = seq_tr - mt_pulse_dur - mt_pulse_int - final_spoiler - ex_dur; % post-readout dead time
readout_par_container{1} = b1_amp_ex;
readout_par_container{2} = N_ex;
readout_par_container{3} = seq_tr;
readout_par_container{4} = ex_dur;
readout_par_container{5} = read_dur;
readout_par_container{6} = seq_dead_time;
readout_par_container{7} = N_shots;
readout_par_container{8} = ex_rf_pulse_profile;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pools_copy = pools;

if length(freq_offsets)>length(pools)
    colors = distinguishable_colors(length(freq_offsets)+N_dummy_offsets);
else
    colors = distinguishable_colors(length(pools));
end

legend_pools = cell(length(pools), 1);
for ii = 1:length(pools)
    legend_pools{ii} = pools{ii}{1};
end

if strcmp(b_mz_evolution, 'false') == 1
    [freq_offsets, spectrum] = cest_slave_sim_ss_fast(spacing_for_mz, pools, freq_offsets, sampling_strategy, ...
        main_field_strength, b1_amp_mt, mt_pulse_dur, mt_pulse_int, final_spoiler,N_seq_repeatitions, N_dummy_offsets, ...
        t1_recovery, readout_par_container, other_seq_par_container);
    
    h1 = figure;
    set(h1, 'name', sprintf(strcat('cest_sequence_', cest_sequence_type, '_sampling_pattern_', sampling_strategy)), 'numbertitle', 'off')
    
    % CEST spectrum
    subplot(121)
    plot(freq_offsets, spectrum, '.-', 'MarkerSize', 30), xlim([-10 10]), ylim([0 1.05])
    set(gca, 'Xdir', 'reverse')
    hold on, grid on
    xlabel('offset [ppm]')
    ylabel('S/S0')
    title('CEST spectrum')
    
    % CEST spectrum decomposition
    subplot(122)
    Za(:,length(pools)+1)=spectrum;
    %each single pool except for free H2O with conc and exch rate set to 0
    for ii=2:length(pools)
        pools=pools_copy;
        pools{ii}{5}=0;
        pools{ii}{6}=0;
        [~, spectrum] = cest_slave_sim_ss_fast(spacing_for_mz, pools, freq_offsets, sampling_strategy, ...
            main_field_strength, b1_amp_mt, mt_pulse_dur, mt_pulse_int, final_spoiler,N_seq_repeatitions, N_dummy_offsets, ...
            t1_recovery, readout_par_container, other_seq_par_container);
        Za(:,ii) = spectrum;
    end
    
    % free H2O
    pools=pools_copy;
    for ii=2:length(pools)
        pools{ii}{5}=0;
        pools{ii}{6}=0;
    end
    [~, spectrum] = cest_slave_sim_ss_fast(spacing_for_mz, pools, freq_offsets, sampling_strategy, ...
        main_field_strength, b1_amp_mt, mt_pulse_dur, mt_pulse_int, final_spoiler,N_seq_repeatitions, N_dummy_offsets, ...
        t1_recovery, readout_par_container, other_seq_par_container);
    Za(:,1) = spectrum;
    
    % pure signal for each pool
    pureZa(:,1)=1-Za(:,1);
    
    for ii=2:length(pools)
        pureZa(:,ii)=Za(:,ii)-Za(:,end);
    end
    
    for ii=1:length(pools)
        plot(freq_offsets, pureZa(:,ii),'Color',colors(ii,:),'LineWidth',4), xlim([-10 10]), ylim([0 1.05])
        hold on, grid on
    end
    set(gca, 'Xdir', 'reverse')   
    xlabel('offset [ppm]')
    ylabel('S/S0')
    title('CEST spectrum: pool decomposition')
    legend(legend_pools)
else
    [freq_offsets, spectrum, offset_order, time_evolution, spectrum_evolution] = cest_slave_sim_ss_evolution(spacing_for_mz, pools, freq_offsets, sampling_strategy, ...
        main_field_strength, b1_amp_mt, mt_pulse_dur, mt_pulse_int, final_spoiler,N_seq_repeatitions, N_dummy_offsets, ...
        t1_recovery, readout_par_container, other_seq_par_container);
    
    h1=figure;
    set(h1, 'name', sprintf(strcat('cest_sequence_', cest_sequence_type, '_sampling_pattern_', sampling_strategy)), 'numbertitle', 'off')
    
    % CEST spectrum
    subplot(121)
    plot(freq_offsets, spectrum, '.-', 'MarkerSize', 30), xlim([-10 10]), ylim([0 1.05])
    set(gca, 'Xdir', 'reverse')
    hold on, grid on
    xlabel('offset [ppm]')
    ylabel('S/S0')
    title('CEST spectrum')

    % CEST spectrum, Mz evolution
    subplot(122)
    for ii=1:length(freq_offsets)
        plot(time_evolution(ii,:), spectrum_evolution(ii,:),'Color',colors(ii,:))
        hold on, grid on
        
        % freq offsets labels
        text(time_evolution(ii,round(length(time_evolution(ii,:))/5)), 1.05, num2str(offset_order(ii)),'Color',colors(ii,:))
    end
    ylim([0 1.05])
    xlabel('Time [s]')
    ylabel('S/S0')
    title('CEST spectrum: Mz evolution')
    
    if length(freq_offsets)<=20
        legend_offsets=cell(1,length(freq_offsets));
        for ii=1:length(freq_offsets)
            legend_offsets{ii}=sprintf('%0.2f ppm (offset# %d)',freq_offsets(ii),offset_order(ii));
        end
        legend(legend_offsets, 'location', 'southeast')
    end
end

end






