function [freq_offsets, spectrum] = cest_slave_sim_ss_fast(~, pools, freq_offsets, sampling_strategy, main_field_strength, ...
    b1_amp_mt, mt_pulse_dur, mt_pulse_int, final_spoiler,N_seq_repeatitions, N_dummy_offsets, t1_recovery, ...
    readout_par_container, other_seq_par_container)

simulation_timer = tic;

% setting up freq offsets based on the sampling pattern and adding (optional) dummies
freq_offsets = rearrange_offsets(sampling_strategy, N_dummy_offsets, freq_offsets);

% setting up all parameters
k0_time = other_seq_par_container{1};
ds_time = other_seq_par_container{2};
mt_rf_pulse_profile = other_seq_par_container{3};
b1_amp_ex = readout_par_container{1};
N_ex = readout_par_container{2};
seq_tr = readout_par_container{3};
ex_dur = readout_par_container{4};
read_dur = readout_par_container{5};
seq_dead_time = readout_par_container{6};
N_shots = readout_par_container{7};
ex_rf_pulse_profile = readout_par_container{8};

N_points_mt = length(mt_rf_pulse_profile);
sampling_rate_mt = mt_pulse_dur/N_points_mt;
N_points_ex=length(ex_rf_pulse_profile);
sampling_rate_ex=ex_dur/N_points_ex;

%% saturation
pulse_cest_mt(1:N_points_mt,1)=mt_rf_pulse_profile./max(mt_rf_pulse_profile)*b1_amp_mt*42.576;
pulse_cest_mt(1:N_points_mt,2)=sampling_rate_mt;
%% interval + spoiler
pulse_cest_int(1,1)=0;
pulse_cest_int(1,2)=mt_pulse_int+final_spoiler;
%% readout
pulse_readout(1:N_points_ex,1)=ex_rf_pulse_profile./max(ex_rf_pulse_profile)*b1_amp_ex*42.576;
pulse_readout(1:N_points_ex,2)=sampling_rate_ex;
%% dead time after readout
pulse_postreadout_dead_time(1,1)=0;
pulse_postreadout_dead_time(1,2)=seq_dead_time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name=cell(1,length(pools));
t1=zeros(1,length(pools));
t2=zeros(1,length(pools));
w_pool_ppm=zeros(1,length(pools));
Con=zeros(1,length(pools));
Rate=zeros(1,length(pools));
lineshape=cell(1,length(pools));
C=zeros(1,length(pools));
Ca=zeros(1,length(pools)-1);
k1=zeros(1,length(pools));
k2=zeros(1,length(pools));

for ii=1:length(pools)
    name{ii}=pools{ii}{1};
    t1(ii)=pools{ii}{2};
    t2(ii)=pools{ii}{3};
    w_pool_ppm(ii)=pools{ii}{4};
    Con(ii)=pools{ii}{5};
    Rate(ii)=pools{ii}{6};
    lineshape{ii}=pools{ii}{7};
end

Con(:)=Con(:)./Con(1);
C(2:end)=Rate(2:end);
Ca=Con.*C;
C(1)=sum(Ca);
k1=1./t1+C;
k2=1./t2+C;
C=C(2:end);
Ca=Ca(2:end);

% w should be in rad/s
ww=freq_offsets*main_field_strength*42.576*2*pi;
w_pool=w_pool_ppm*main_field_strength*42.576*2*pi;

spectrum=zeros(length(freq_offsets),1);
y=[zeros(1,length(pools)*2) Con]';
b=[zeros(1,length(pools)*2) Con./t1]';

for offset_No=1:length(ww)
    
    w=ww(offset_No);
    
    %% saturation
    expmA_cest_dur = zeros(length(pulse_cest_mt),length(pools)*3,length(pools)*3);
    Ainvbm_cest_dur = zeros(length(pulse_cest_mt), length(pools)*3);
    
    for m=1:size(pulse_cest_mt,1)
        w1=2.0*pi*pulse_cest_mt(m,1);
        add_lineshape=zeros(1,length(pools));
        for ii=1:length(pools)
            if isempty(lineshape{ii}) || strcmp(lineshape{ii},'lorentzian')==1
                add_lineshape(ii)=0;
            elseif strcmp(lineshape{ii},'gaussian')==1
                gauss=t2(ii)/sqrt(2*pi)*exp( -2.0 * ( pi * (w-w_pool(ii)) * t2(ii))^2.0 );
                add_lineshape(ii)=pi*w1*w1*gauss;
            elseif strcmp(lineshape{ii},'super_lorentzian')==1
                th=linspace(0,pi/2,1000);
                Y=@(th) sin(th).*sqrt(2/pi)*t2(ii)./abs(3*cos(th).*cos(th)-1).*exp(-2*(((w-w_pool(ii)))*t2(ii)./abs(3*cos(th).*cos(th)-1)).^2);
                s_lor=integral(Y,0,pi/2,'AbsTol',1e-9);
                add_lineshape(ii)=pi*w1*w1*s_lor;
            end
        end
        
        % diagonal
        A=zeros(length(pools)*3,length(pools)*3);
        A=A+diag([-k2 -k2 -(k1+add_lineshape)]);
        
        % off diagonals
        for ii=1:length(pools)-1
            upperC=repmat([C(ii) zeros(1,length(pools)-1)],[1,3]);
            upperC=upperC(1:end-ii);
            A=A+diag(upperC,ii);
            
            lowerC=repmat([Ca(ii) zeros(1,length(pools)-1)],[1,3]);
            lowerC=lowerC(1:end-ii);
            A=A+diag(lowerC,-ii);
        end
        A=A+diag([-(w_pool-w) repmat(-w1,[1,length(pools)])],length(pools));
        A=A+diag([(w_pool-w) repmat(w1,[1,length(pools)])],-length(pools));
        
        expmA_cest_dur(m,:,:)= expm(A*pulse_cest_mt(m,2));
        Ainvbm_cest_dur(m,:) = A\b;
    end
    
    %% interval + spoiler
    expmA_cest_int = zeros(length(pulse_cest_int),length(pools)*3,length(pools)*3);
    Ainvbm_cest_int = zeros(length(pulse_cest_int), length(pools)*3);
    
    for m=1:size(pulse_cest_int,1)
        w1=2.0*pi*pulse_cest_int(m,1);
        add_lineshape=zeros(1,length(pools));
        
        % diagonal
        A=zeros(length(pools)*3,length(pools)*3);
        A=A+diag([-k2 -k2 -(k1+add_lineshape)]);
        
        % off diagonals
        for ii=1:length(pools)-1
            upperC=repmat([C(ii) zeros(1,length(pools)-1)],[1,3]);
            upperC=upperC(1:end-ii);
            A=A+diag(upperC,ii);
            
            lowerC=repmat([Ca(ii) zeros(1,length(pools)-1)],[1,3]);
            lowerC=lowerC(1:end-ii);
            A=A+diag(lowerC,-ii);
        end
        A=A+diag([-(w_pool-w) repmat(-w1,[1,length(pools)])],length(pools));
        A=A+diag([(w_pool-w) repmat(w1,[1,length(pools)])],-length(pools));
        
        expmA_cest_int(m,:,:)= expm(A*pulse_cest_int(m,2));
        Ainvbm_cest_int(m,:) = A\b;
    end
       
    %% readout
    expmA_readout = zeros(length(pulse_readout),length(pools)*3,length(pools)*3);
    Ainvbm_readout = zeros(length(pulse_readout), length(pools)*3);
    
    for m=1:size(pulse_readout,1)
        w1=2.0*pi*pulse_readout(m,1);
        add_lineshape=zeros(1,length(pools));
        
        % diagonal
        A=zeros(length(pools)*3,length(pools)*3);
        A=A+diag([-k2 -k2 -(k1+add_lineshape)]);
        
        % off diagonals
        for ii=1:length(pools)-1
            upperC=repmat([C(ii) zeros(1,length(pools)-1)],[1,3]);
            upperC=upperC(1:end-ii);
            A=A+diag(upperC,ii);
            
            lowerC=repmat([Ca(ii) zeros(1,length(pools)-1)],[1,3]);
            lowerC=lowerC(1:end-ii);
            A=A+diag(lowerC,-ii);
        end
        A=A+diag([-(w_pool-0) repmat(-w1,[1,length(pools)])],length(pools));
        A=A+diag([(w_pool-0) repmat(w1,[1,length(pools)])],-length(pools));
        
        expmA_readout(m,:,:)= expm(A*pulse_readout(m,2));
        Ainvbm_readout(m,:) = A\b;
    end
    
    %% post-readout dead time
    expmA_deadtime = zeros(length(pulse_postreadout_dead_time),length(pools)*3,length(pools)*3);
    Ainvbm_deadtime = zeros(length(pulse_postreadout_dead_time), length(pools)*3);
    
    for m=1:size(pulse_postreadout_dead_time,1)
        w1=2.0*pi*pulse_postreadout_dead_time(m,1);
        add_lineshape=zeros(1,length(pools));
        
        % diagonal
        A=zeros(length(pools)*3,length(pools)*3);
        A=A+diag([-k2 -k2 -(k1+add_lineshape)]);
        
        % off diagonals
        for ii=1:length(pools)-1
            upperC=repmat([C(ii) zeros(1,length(pools)-1)],[1,3]);
            upperC=upperC(1:end-ii);
            A=A+diag(upperC,ii);
            
            lowerC=repmat([Ca(ii) zeros(1,length(pools)-1)],[1,3]);
            lowerC=lowerC(1:end-ii);
            A=A+diag(lowerC,-ii);
        end
        A=A+diag([-(w_pool-w) repmat(-w1,[1,length(pools)])],length(pools));
        A=A+diag([(w_pool-w) repmat(w1,[1,length(pools)])],-length(pools));
        
        expmA_deadtime(m,:,:)= expm(A*pulse_postreadout_dead_time(m,2));
        Ainvbm_deadtime(m,:) = A\b;
    end
    
    %% inter volume T1 recovery
    w1=0;
    % diagonal
    A=zeros(length(pools)*3,length(pools)*3);
    A=A+diag([-k2 -k2 -(k1+add_lineshape)]);
    
    % off diagonals
    for ii=1:length(pools)-1
        upperC=repmat([C(ii) zeros(1,length(pools)-1)],[1,3]);
        upperC=upperC(1:end-ii);
        A=A+diag(upperC,ii);
        
        lowerC=repmat([Ca(ii) zeros(1,length(pools)-1)],[1,3]);
        lowerC=lowerC(1:end-ii);
        A=A+diag(lowerC,-ii);
    end
    A=A+diag([-(w_pool-w) repmat(-w1,[1,length(pools)])],length(pools));
    A=A+diag([(w_pool-w) repmat(w1,[1,length(pools)])],-length(pools));
    
    expmA_rec = expm(A*t1_recovery);
    Ainvbm_rec= A\b;
    
    for seq_repeat_No=1:N_seq_repeatitions
        %% saturation
        for m=1:size(pulse_cest_mt,1)
            y = reshape(expmA_cest_dur(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_cest_dur(m,:)')-Ainvbm_cest_dur(m,:)';
        end
        
        %% interval + spoiler
        for m=1:size(pulse_cest_int,1)
            y(1:length(pools)*2)=0; % spoiling
            y = reshape(expmA_cest_int(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_cest_int(m,:)')-Ainvbm_cest_int(m,:)';
        end
        
        y(1:length(pools)*2)=0; % spoiling        
        if seq_repeat_No == ceil(k0_time/seq_tr) % k-space center
            spectrum(offset_No) = abs(y(length(pools)*2+1));
        end
        
        %% readout
        for ex_No=1:N_ex
            for m=1:size(pulse_readout,1)
                y = reshape(expmA_readout(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_readout(m,:)')-Ainvbm_readout(m,:)';
            end
            y(1:length(pools)*2)=0; % spoiling
        end
        y(1:length(pools)*2)=0; % spoiling
        
        %% post-readout dead time
        for m=1:size(pulse_postreadout_dead_time,1)
            y = reshape(expmA_deadtime(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_deadtime(m,:)')-Ainvbm_deadtime(m,:)';
        end
    end
    
    %% intervolume T1 recovery
    if t1_recovery>0
        y = expmA_rec*(y+Ainvbm_rec)-Ainvbm_rec;
    end
end

% removing dummies and correcting freq offsets
[freq_offsets, spectrum] = undo_rearrange_offsets(N_dummy_offsets, freq_offsets, spectrum);

elapsed_time = toc(simulation_timer);
fprintf('Elapsed time is %d minutes and %f seconds\n', floor(elapsed_time/60), rem(elapsed_time, 60))
end