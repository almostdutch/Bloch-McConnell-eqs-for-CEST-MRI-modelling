function [freq_offsets, spectrum] = undo_rearrange_offsets(N_dummy_offsets, freq_offsets, spectrum)

spectrum=spectrum(1+N_dummy_offsets:end); % remove dummies
freq_offsets=freq_offsets(1+N_dummy_offsets:end); % remove dummies
[freq_offsets, offset_order]=sort(freq_offsets);
spectrum=spectrum(offset_order);

freq_offsets = freq_offsets(:);
spectrum = spectrum(:);
end