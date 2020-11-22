function [freq_offsets] = rearrange_offsets(sampling_strategy, N_dummy_offsets, freq_offsets)

if strcmp(sampling_strategy, 'jump')
    pos=find(freq_offsets>0);
    neg=find(freq_offsets<0);
    new_offset=zeros(size(freq_offsets));
    new_offset(1)=0;
    new_offset([2:2:length(freq_offsets)])=flipdim(freq_offsets(neg(:)),2);
    new_offset([3:2:length(freq_offsets)])=freq_offsets(pos(:));
    freq_offsets=[zeros(1,N_dummy_offsets) new_offset(:)']; % add dumies
elseif strcmp(sampling_strategy, 'linear')
    freq_offsets=[repmat(freq_offsets(1),[1 N_dummy_offsets]) freq_offsets(:)']; % add dummies
elseif strcmp(sampling_strategy,'random')
    freq_offsets=freq_offsets(randperm(length(freq_offsets)));
    freq_offsets=[repmat(freq_offsets(1),[1 N_dummy_offsets]) freq_offsets(:)']; % add dumies
end

end