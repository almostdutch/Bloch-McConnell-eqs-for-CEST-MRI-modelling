function [freq_offsets, spectrum, offset_order, time_evolution, spectrum_evolution] = undo_rearrange_offsets_evolution(N_dummy_offsets, freq_offsets, spectrum, time_evolution, spectrum_evolution)

spectrum=spectrum(1+N_dummy_offsets:end); % remove dummies
freq_offsets=freq_offsets(1+N_dummy_offsets:end); % remove dummies
time_evolution = time_evolution(1+N_dummy_offsets:end,:); % remove dummies
spectrum_evolution = spectrum_evolution(1+N_dummy_offsets:end,:); % remove dummies

[freq_offsets, offset_order]=sort(freq_offsets);
spectrum=spectrum(offset_order);
time_evolution = time_evolution(offset_order,:);
spectrum_evolution = spectrum_evolution(offset_order,:);

freq_offsets = freq_offsets(:);
spectrum = spectrum(:);

end