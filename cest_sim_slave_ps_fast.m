function varargout=cest_sim_slave_ps_fast(spacing, pools, offset_w, pulse_method, cwpe_approximation, sampling, fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par)
%function varargout=cest_sim_slave_ss(spacing, pools, offset_w, pulse_method, cwpe_approximation, sampling, fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par)

if strcmp(pulse_method,'ps')~=1
    sprintf('wrong_cest_sim_file')
    stop
end

simulation_start=tic;

if strcmp(sampling,'jump')
    pos=find(offset_w>0);
    neg=find(offset_w<0);
    new_offset=zeros(size(offset_w));
    new_offset(1)=0;
    new_offset([2:2:length(offset_w)])=flipdim(offset_w(neg(:)),2);
    new_offset([3:2:length(offset_w)])=offset_w(pos(:));
    offset=[zeros(1,dummy) new_offset(:)'];
elseif strcmp(sampling,'linear')
    offset=[repmat(offset_w(1),[1 dummy]) offset_w(:)'];
elseif strcmp(sampling,'random')
    offset=offset_w(randperm(length(offset_w)));
    offset=[repmat(offset(1),[1 dummy]) offset(:)'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CEST saturation
am_sg_100_100_0 = [0, 413, 859, 1337, 1847, 2389, 2962, 3568, 4204, 4870,...
    5565, 6288, 7038, 7814, 8613, 9435, 10277, 11138, 12015, 12906,...
    13809, 14721, 15640, 16563, 17487, 18409, 19327, 20237, 21136, 22023,...
    22892, 23742, 24570, 25372, 26146, 26889, 27598, 28271, 28905, 29497,...
    30046, 30549, 31005, 31411, 31767, 32070, 32320, 32515, 32655, 32739,...
    32767, 32739, 32655, 32515, 32320, 32070, 31767, 31411, 31005, 30549,...
    30046, 29497, 28905, 28271, 27598, 26889, 26146, 25372, 24570, 23742,...
    22892, 22023, 21136, 20237, 19327, 18409, 17487, 16563, 15640, 14721,...
    13809, 12906, 12015, 11138, 10277, 9435, 8613, 7814, 7038, 6288,...
    5565, 4870, 4204, 3568, 2962, 2389, 1847, 1337, 859, 413, 0 ];

saturation_RF=am_sg_100_100_0;
saturation_RF=ones(size(am_sg_100_100_0)); %% equivalent of a perfect rectangular
%pulse

%% up/down sampling the saturation prepulse
% saturation_RF=interp(saturation_RF,pulse_sampling_factor);
% saturation_RF=saturation_RF(1:end-pulse_sampling_factor+1);

%% excitation
am_sg_400_150_125 =[0, 47, 96, 147, 200, 253, 305, 355, 403, 448,...
    487, 521, 548, 567, 577, 578, 569, 549, 518, 476,...
    422, 358, 283, 197, 103, 0, -110, -225, -343, -464,...
    -584, -703, -816, -923, -1021, -1108, -1181, -1237, -1277, -1296,...
    -1294, -1270, -1222, -1150, -1053, -933, -788, -621, -432, -225,...
    0, 239, 488, 744, 1003, 1262, 1514, 1756, 1983, 2190,...
    2372, 2524, 2643, 2724, 2763, 2757, 2703, 2600, 2446, 2240,...
    1983, 1676, 1320, 920, 478, 0, -509, -1042, -1593, -2152,...
    -2712, -3263, -3795, -4298, -4763, -5178, -5535, -5821, -6029, -6150,...
    -6173, -6093, -5902, -5595, -5168, -4617, -3941, -3139, -2213, -1165,...
    0, 1276, 2657, 4132, 5692, 7326, 9021, 10762, 12536, 14327,...
    16119, 17896, 19641, 21338, 22970, 24522, 25977, 27323, 28544, 29629,...
    30567, 31348, 31964, 32409, 32677, 32767, 32677, 32409, 31964, 31348,...
    30567, 29629, 28544, 27323, 25977, 24522, 22970, 21338, 19641, 17896,...
    16119, 14327, 12536, 10762, 9021, 7326, 5692, 4132, 2657, 1276,...
    0, -1165, -2213, -3139, -3941, -4617, -5168, -5595, -5902, -6093,...
    -6173, -6150, -6029, -5821, -5535, -5178, -4763, -4298, -3795, -3263,...
    -2712, -2152, -1593, -1042, -509, 0, 478, 920, 1320, 1676,...
    1983, 2240, 2446, 2600, 2703, 2757, 2763, 2724, 2643, 2524,...
    2372, 2190, 1983, 1756, 1514, 1262, 1003, 744, 488, 239,...
    0];

% am_sg_150_100_167=[66, -138, -359, -595, -845, -1107, -1378, -1658, -1943, -2230,...
%     -2517, -2800, -3076, -3341, -3592, -3824, -4034, -4217, -4370, -4488,...
%     -4567, -4602, -4591, -4530, -4414, -4240, -4007, -3710, -3348, -2918,...
%     -2421, -1853, -1216, -510, 266, 1110, 2019, 2991, 4024, 5114,...
%     6256, 7446, 8678, 9947, 11247, 12571, 13911, 15262, 16614, 17960,...
%     19293, 20604, 21885, 23129, 24326, 25471, 26554, 27570, 28511, 29371,...
%     30144, 30826, 31410, 31894, 32274, 32547, 32712, 32767, 32712, 32547,...
%     32274, 31894, 31410, 30826, 30144, 29371, 28511, 27570, 26554, 25471,...
%     24326, 23129, 21885, 20604, 19293, 17960, 16614, 15262, 13911, 12571,...
%     11247, 9947, 8678, 7446, 6256, 5114, 4024, 2991, 2019, 1110,...
%     266];
excitation_RF=am_sg_400_150_125;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sat_amp=B1amp;
ex_amp=readout_par{1};
NofEx = readout_par{2};
TR1=readout_par{3};
ex_dur=readout_par{4};
TR2=readout_par{5};
after_TR_dead_time=readout_par{6};
NofShots=readout_par{7};

sampling_rate_sat=pulse_dur/length(saturation_RF);
Nsat=length(saturation_RF);
Nint=round(pulse_int/sampling_rate_sat);
Nsp=round(spoiler/sampling_rate_sat);
sampling_rate_ex=ex_dur/length(excitation_RF);
Nex=length(excitation_RF);
Nrd=round((TR2-ex_dur)/sampling_rate_ex);
Ndeadtime = round(after_TR_dead_time/sampling_rate_sat);

%% saturation
if strcmp(cwpe_approximation,'false')==1
    pulse_cest_dur(1:Nsat,1)=saturation_RF./max(saturation_RF)*sat_amp*42.576;
    pulse_cest_dur(1:Nsat,2)=sampling_rate_sat;
elseif strcmp(cwpe_approximation,'true')==1
    pulse_cest_dur(1,1)=sat_amp*42.576;
    pulse_cest_dur(1,2)=pulse_dur;
end
%% interval
pulse_cest_int(1,1)=0;
pulse_cest_int(1,2)=pulse_int;
%% final spoiler
pulse_sp(1,1)=0;
pulse_sp(1,2)=spoiler;
%% readout
pulse_readout(1:Nex,1)=excitation_RF./max(excitation_RF)*ex_amp*42.576;
pulse_readout(1:Nex,2)=sampling_rate_ex;
%% post-excitation dead time
pulse_postexcitation_dead_time(1,1)=0;
pulse_postexcitation_dead_time(1,2)=TR2-ex_dur;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

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
        if strcmp(name{ii},'MT')==1
        sl_n=ii;
    end

    t1(ii)=pools{ii}{2};
    t2(ii)=pools{ii}{3};
    w_pool_ppm(ii)=pools{ii}{4};
    Con(ii)=pools{ii}{5};
    Rate(ii)=pools{ii}{6};
    lineshape{ii}=pools{ii}{7};
end

Con(:)=Con(:)./Con(1); % normalized conc
C(2:end)=Rate(2:end);
Ca=Con.*C;
C(1)=sum(Ca);
k1=1./t1+C;
k2=1./t2+C;
C=C(2:end);
Ca=Ca(2:end);

% w should be in rad/s
ww=offset.*fs.*42.576.*2.*pi;
w_pool=w_pool_ppm*fs*42.576*2*pi;
delta=2*fs*42.576*2*pi;
% [~, n_sl1]=min(abs(ww-(w_pool(sl_n)-delta)));
% [~, n_sl2]=min(abs(ww-(w_pool(sl_n)+delta)));
Za=zeros(length(offset),1);
y=[zeros(1,length(pools)*2) Con]';
b=[zeros(1,length(pools)*2) Con./t1]';

for n=1:length(ww);
    
    w=ww(n);
    
    %% saturation
    expmA_cest_dur = zeros(length(pulse_cest_dur),length(pools)*3,length(pools)*3);
    Ainvbm_cest_dur = zeros(length(pulse_cest_dur), length(pools)*3);
    
    for m=1:size(pulse_cest_dur,1);
        w1=2.0*pi*pulse_cest_dur(m,1);
        add_lineshape=zeros(1,length(pools));
        
%                 for ii=1:length(pools)
%                     if isempty(lineshape{ii}) || strcmp(lineshape{ii},'lorentzian')==1
%                         add_lineshape(ii)=0;
%                     elseif strcmp(lineshape{ii},'gaussian')==1
%                         gauss=t2(ii)/sqrt(2*pi)*exp( -2.0 * ( pi * (w-w_pool(ii)) * t2(ii))^2.0 );
%                         add_lineshape(ii)=pi*w1*w1*gauss;
%                     elseif strcmp(lineshape{ii},'super_lorentzian')==1
%                         th=linspace(0,pi/2,1000);
%                         Y=@(th) sin(th).*sqrt(2/pi)*t2(ii)./abs(3*cos(th).*cos(th)-1).*exp(-2*(((w-w_pool(ii)))*t2(ii)./abs(3*cos(th).*cos(th)-1)).^2);
%                         s_lor=integral(Y,0,pi/2,'AbsTol',1e-9);
%                         add_lineshape(ii)=pi*w1*w1*s_lor;
%                     end
%                 end

%         for ii=1:length(pools)
%             if strcmp(lineshape{ii},'super_lorentzian')==1
%                 Y=@(th) sin(th).*sqrt(2/pi)*t2(ii)./abs(3*cos(th).*cos(th)-1).*exp(-2*(((w-w_pool(ii)))*t2(ii)./abs(3*cos(th).*cos(th)-1)).^2);
%                 s_lor=integral(Y,0,pi/2,'AbsTol',1e-9);
%                 
%                 if n==n_sl1
%                     sl_th=s_lor;
%                 end
%                 
%                 if n>n_sl1 && n<n_sl2
%                     s_lor=sl_th;
%                 end
%                 add_lineshape(ii)=pi*w1*w1*s_lor;
%             end
%         end

        
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
        
        expmA_cest_dur(m,:,:)= expm(A*pulse_cest_dur(m,2));
        Ainvbm_cest_dur(m,:) = A\b;
    end
    
    %% interval
    expmA_cest_int = zeros(length(pulse_cest_int),length(pools)*3,length(pools)*3);
    Ainvbm_cest_int = zeros(length(pulse_cest_int), length(pools)*3);
    
    for m=1:size(pulse_cest_int,1);
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
    
    %% final spoiler
    expmA_sp = zeros(length(pulse_sp),length(pools)*3,length(pools)*3);
    Ainvbm_sp = zeros(length(pulse_sp), length(pools)*3);
    
    for m=1:size(pulse_sp,1);
        w1=2.0*pi*pulse_sp(m,1);
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
        
        expmA_sp(m,:,:)= expm(A*pulse_sp(m,2));
        Ainvbm_sp(m,:) = A\b;
    end
    
    %% readout
    expmA_readout = zeros(length(pulse_readout),length(pools)*3,length(pools)*3);
    Ainvbm_readout = zeros(length(pulse_readout), length(pools)*3);
    
    for m=1:size(pulse_readout,1);
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
    expmA_deadtime = zeros(length(pulse_postexcitation_dead_time),length(pools)*3,length(pools)*3);
    Ainvbm_deadtime = zeros(length(pulse_postexcitation_dead_time), length(pools)*3);
    
    for m=1:size(pulse_postexcitation_dead_time,1);
        w1=2.0*pi*pulse_postexcitation_dead_time(m,1);
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
        
        expmA_deadtime(m,:,:)= expm(A*pulse_postexcitation_dead_time(m,2));
        Ainvbm_deadtime(m,:) = A\b;
    end
    
    %% inter shot T1 recovery
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
    
    expmA_rec = expm(A*recovery);
    Ainvbm_rec= A\b;
    
    for N_shots=1:NofShots
        %% saturation
        for N_pulse_repeat=1:pulse_repeat
            
            for m=1:size(pulse_cest_dur,1)
                y = reshape(expmA_cest_dur(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_cest_dur(m,:)')-Ainvbm_cest_dur(m,:)';
            end
            
            %% interval
            if N_pulse_repeat==pulse_repeat
                continue
            end
            
            if pulse_int>0
                for m=1:size(pulse_cest_int,1)
                    y(1:length(pools)*2)=0; % spoiling
                    y = reshape(expmA_cest_int(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_cest_int(m,:)')-Ainvbm_cest_int(m,:)';
                end
            end
            
        end
        
        %% final spoiler
        if pulse_sp>0
            for m=1:size(pulse_sp,1)
                y(1:length(pools)*2)=0; % spoiling
                y = reshape(expmA_sp(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_sp(m,:)')-Ainvbm_sp(m,:)';
            end
        end
        
        Za(n) = abs(y(length(pools)*2+1))*100;
        
        %% readout
        for N_ex=1:NofEx
            for m=1:size(pulse_readout,1)
                y = reshape(expmA_readout(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_readout(m,:)')-Ainvbm_readout(m,:)';
            end
            y(1:length(pools)*2)=0; % spoiling
            
            %% post-excitation dead time
            for m=1:size(pulse_postexcitation_dead_time,1)
                y = reshape(expmA_deadtime(m,:,:),[length(pools)*3, length(pools)*3])*(y+Ainvbm_deadtime(m,:)')-Ainvbm_deadtime(m,:)';
            end
        end
        
        %% intershot T1 recovery
        if recovery>0
            y = expmA_rec*(y+Ainvbm_rec)-Ainvbm_rec;
        end
    end
end


if strcmp(sampling,'jump')
    Za=Za(1+dummy:end);
    offset=offset(1+dummy:end);
    [offset, offset_order]=sort(offset);
    Za=Za(offset_order);
elseif strcmp(sampling,'linear')
    Za=Za(1+dummy:end);
    offset=offset(1+dummy:end);
elseif strcmp(sampling,'random')
    Za=Za(1+dummy:end);
    offset=offset(1+dummy:end);
    save offset_random
    [offset, offset_order]=sort(offset);
    Za=Za(offset_order);
end

% [~, n_sl1]=min(abs(offset-(w_pool_ppm(sl_n)-8)));
% [~, n_sl2]=min(abs(offset-(w_pool_ppm(sl_n)+2)));
% 
% Za(n_sl1:n_sl2) = smooth(Za(n_sl1:n_sl2),0.2,'rloess');
varargout{1}=Za;
varargout{2}=offset;


ElapsedScript=toc(simulation_start);
fprintf('Elapsed time is %d minutes and %f seconds\n', floor(ElapsedScript/60),rem(ElapsedScript,60))