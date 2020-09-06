function cest_master

%close all
%amp2cwpe_ps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
tissue='GM';
sequence='ps'; % ss - steady state (also interleaved), ps - pseudo-steady state (also long)
Mz_evolution=0;
single_pool_decomposition=0;

if strcmp(tissue,'WM')==1
    %WM
    pools={...
        {'H2O',1.1,   0.040, 0, 2*55.6*1000, 0,'lorentzian'}...
        {'APT',1.0, 10e-03, 3.5, 166.8, 50,'lorentzian'}...
        {'MT',1.0, 10e-06, -5, 12232, 50,'super_lorentzian'}...
        {'NOE',1.0, 0.3e-03, -3,5, 6672, 10,'lorentzian'}...
        };
elseif strcmp(tissue,'GM')==1
    %GM
    pools={...
        {'H2O',1.9,   0.055, 0, 2*55.6*1000, 0,'lorentzian'}...
        {'APT',1.0, 10e-03, 3.5, 166.8, 300,'lorentzian'}...
        {'MT',1.0, 10e-06, -2.4, 14456, 50,'lorentzian'}...
        {'NOE',1.0, 0.3e-03, -3.5, 3336, 10,'lorentzian'}...
        };
end
% 
% pools={...
%     {'H2O',4,   0.50, 0, 2*55.6*1000, 0,'lorentzian'}...
%     {'GAG',1.0, 10e-03, 2, 2*55.6*1000*0.27/100, 1000,'lorentzian'}...
%     };

% pools={...
%     {'H2O',1.8,   0.040, 0, 2*55.6*1000, 0,'lorentzian'}...
%     {'MT',1.0, 10e-06, 0, 5000, 10,'super_lorentzian'}...
%     {'cr',1.0, 10e-03, 2, 850, 50,'lorentzian'}...
%     {'biomat',1.0, 10e-03, 5, 1000, 10,'lorentzian'}...
%     {'cells',1.0, 10e-03, 30, 500, 10,'lorentzian'}...
%     };

offset_w=[-1.5:0.1:1.5]; % offsets in ppm, the last one is for normalization
%% pulse sequence parameters
if strcmp(sequence,'ss')==1
    pulse_method='ss';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this parameters can be changed
    cwpe_approximation='false'; % true or false
    sampling='linear'; % jump, linear or random (only linear for now)
    fs=7.0;
    B1amp= 5; % added to cest_sim_slace_ss_fast.m for cwpe /1.2987 (found this by playing with simulations)
    pulse_dur=5/1000;
    pulse_int=0;
    spoiler=0.5/1000;
    TR1=65/1000; % TR
    k0_time=8.8; % in s
    ds_time=16.9; % in s
    pulse_repeat=ceil(ds_time/TR1);
    dummy=0;
    recovery=10; % in s, inter volume T1 recovery
    spacing=ds_time+recovery; % in s, just spacing out ss progress
    %% readout
    % EPI or TFE
    B1amp_ex=1.521; % in uT
    NofEx = 1;	% excitations.
    NofShots=1;
    ex_dur=0.2240/1000; % in s
    TR2=8.3/1000; % readout duration
    pulse_sampling_factor=10; % determines sampling rate
elseif strcmp(sequence,'ps')==1
    pulse_method='ps';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this parameters can be changed
    cwpe_approximation='false'; % false or true
    sampling='linear'; % jump, linear or random
    fs=7;
    B1amp=3.7;%0.8978;  rect 2ut (40, dc100%) - 2.7818 (20, dc50%) - 2.0006 (1, dc100%); sinc 2ut amp (40, dc100%) - 2.7818ut amp (20, dc50%) - 1.2669 avg (1, dc100%)
    pulse_dur=5/1000;
    pulse_int=0/1000;
    spoiler=0.4/1000;
    TR1=0.000; % not used
    pulse_repeat=1;
    dummy=0;
    recovery=15; % in s, inter shot T1 recovery
    pulse_sampling_factor=1; % determines sampling rate
    %% readout
    % EPI or TFE
    B1amp_ex=5; % in uT
    NofEx = 10;	% excitations.
    NofShots=1;
    ex_dur=0.2/1000; % in s
    TR2=2.5/1000; % TR of readout
    spacing=(pulse_dur*pulse_repeat+pulse_int*(pulse_repeat-1)+spoiler+(ex_dur+TR2)*NofEx+recovery)*NofShots; % in s, just spacing out ss progress
end
% figure
%% additional parameters
other_par=[];
if exist('k0_time','var')==1 && exist('ds_time','var')==1
    other_par{1}=k0_time;
    other_par{2}=ds_time;
    other_par{3}=Mz_evolution;
end
after_TR_dead_time=TR1-pulse_dur-pulse_int-spoiler-ex_dur; % no used for long
readout_par{1}=B1amp_ex;
readout_par{2}=NofEx;
readout_par{3}=TR1;
readout_par{4}=ex_dur;
readout_par{5}=TR2;
readout_par{6}=after_TR_dead_time;
readout_par{7}=NofShots;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pools_copy=pools;

if length(offset_w)>length(pools)
    colors = distinguishable_colors(length(offset_w)+dummy);
else
    colors = distinguishable_colors(length(pools));
end

legend_pools=cell(1,length(pools));
legend_pools{1}='spectrum';
for ii=1:length(pools)
    legend_pools{ii+1}=pools{ii}{1};
end


if Mz_evolution==1 && single_pool_decomposition==0
    [Za,time_all, sp_all, ~, ~, offset_w,offset_order]=cest_sim_evolution(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
    Za=Za/100;
    sp_all=sp_all/100;
    [Za_norm, ~]=cest_sim_fast(spacing,pools, 500, pulse_method, cwpe_approximation, 'linear', ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
    Za_norm=Za_norm/100;
    Za=Za./Za_norm*100;
    sp_all=sp_all./Za_norm*100;
    
    h1=figure;
    for ii=1:length(offset_w)
        figure(h1)
        plot(time_all(ii,:),sp_all(ii,:),'Color',colors(ii,:))
        hold on
        grid on
        
        %% numbering offsets
        text(time_all(ii,round(length(time_all(ii,:))/5)),105,num2str(offset_order(ii)),'Color',colors(ii,:))
        plot([time_all(ii,round(length(time_all(ii,:))/5)) time_all(ii,round(length(time_all(ii,:))/5))],[max(sp_all(ii,:))+2 105],'Color',colors(ii,:))
    end
    ylim([0 110])
    xlabel('Time, s')
    ylabel('Mo/Mz, %')
    
    set(h1, 'name',sprintf(strcat('mz_evolution','_',pulse_method,'_',sampling)),'numbertitle','off')
    
    if length(offset_w)<=20
        legend_offsets=cell(1,length(offset_w));
        for ii=1:length(offset_w)
            legend_offsets{ii}=sprintf('%0.2f ppm',offset_w(ii));
        end
        figure(h1)
        legend(legend_offsets)
    end
    
    h2=figure;
    Za=Za(:);
    offset_w=offset_w(:);
    figure(h2)
    for ii=1:length(offset_w)
        plot(offset_w(ii),Za(ii),'.','MarkerSize',30,'Color',colors(ii,:))
        hold on,grid on
        text(offset_w(ii),Za(ii)+3,num2str(offset_order(ii)),'Color',colors(ii,:))
        %plot([offset_w(ii) offset_w(ii)],[Za(ii)+2 105],'Color',colors(ii,:))
    end
    set(gca,'Xdir','reverse')
    hold on,grid on
    xlim([-10 10])
    ylim([0 110])
    xlabel('offset')
    ylabel('S/S0')
    set(gca,'fontsize',24)
    set(h2, 'name',sprintf(strcat('cest_color_coding','_',pulse_method,'_',sampling)),'numbertitle','off')
    
    offset_0=round((length(offset_w)-1)/2);
    for tt=1:offset_0
        MTR(tt)=Za(offset_0+1-tt)-Za(offset_0+1+tt);%./Za(offset_0+1-tt);
    end
    MTR_w=offset_w(offset_0+2:end);
    
    h3=figure;
    figure(h3)
    subplot(1,2,1)
    plot(offset_w, Za,'.','MarkerSize',30)
    set(gca,'Xdir','reverse')
    hold on,grid on
    xlim([-10 10])
    ylim([0 110])
    xlabel('offset')
    ylabel('S/S0')
    set(gca,'fontsize',24)
    subplot(1,2,2)
    plot(MTR_w, MTR,'.','MarkerSize',30)
    hold on,grid on
    xlim([0 10])
    ylim([-2 10])
    xlabel('offset')
    ylabel('MTR %')
    set(gca,'fontsize',24)
    set(h3, 'name',sprintf(strcat('cest_mtr','_',pulse_method,'_',sampling)),'numbertitle','off')
end

if Mz_evolution==0 && single_pool_decomposition==0
    [Za, offset_w]=cest_sim_fast(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
    Za=Za/100;
    [Za_norm, ~]=cest_sim_fast(spacing,pools, 500, pulse_method, cwpe_approximation, 'linear', ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
    Za_norm=Za_norm/100;
    
    Za=Za/Za_norm*100;
    Za=Za(:);
    offset_w=offset_w(:);
    
    offset_0=round((length(offset_w)-1)/2);
    for tt=1:offset_0
        MTR(tt)=Za(offset_0+1-tt)-Za(offset_0+1+tt);%./Za(offset_0+1-tt);
    end
    MTR_w=offset_w(offset_0+2:end);
    
    %h3=figure;
    %figure(h3)
    subplot(1,2,1)
    plot(offset_w, Za,'.-','MarkerSize',30)
    set(gca,'Xdir','reverse')
    hold on,grid on
    xlim([-2 2])
    ylim([0 110])
    xlabel('offset')
    ylabel('S/S0')
    set(gca,'fontsize',24)
    subplot(1,2,2)
    plot(MTR_w, MTR,'.','MarkerSize',30)
    hold on,grid on
    xlim([0 10])
    ylim([-2 10])
    xlabel('offset')
    ylabel('MTR %')
    set(gca,'fontsize',24)
    set(h3, 'name',sprintf(strcat('cest_mtr','_',pulse_method,'_',sampling)),'numbertitle','off')
end

% spectrum decomposition
if Mz_evolution==0 && single_pool_decomposition==1
    [Za, offset_w]=cest_sim_fast(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
    Za=Za/100;
    [Za_norm, ~]=cest_sim_fast(spacing,pools, 500, pulse_method, cwpe_approximation, 'linear', ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
    Za_norm=Za_norm/100;
    
    Za=Za/Za_norm*100;
    Za=Za(:);
    offset_w=offset_w(:);
    
    offset_0=round((length(offset_w)-1)/2);
    for tt=1:offset_0
        MTR(tt)=Za(offset_0+1-tt)-Za(offset_0+1+tt);%./Za(offset_0+1-tt);
    end
    MTR_w=offset_w(offset_0+2:end);
    
    h3=figure;
    figure(h3)
    subplot(1,2,1)
    plot(offset_w, Za,'.','MarkerSize',30)
    set(gca,'Xdir','reverse')
    hold on,grid on
    xlim([-5 5])
    ylim([0 110])
    xlabel('offset')
    ylabel('S/S0')
    set(gca,'fontsize',24)
    subplot(1,2,2)
    plot(MTR_w, MTR,'.','MarkerSize',30)
    hold on,grid on
    xlim([0 5])
    ylim([-2 10])
    xlabel('offset')
    ylabel('MTR %')
    set(gca,'fontsize',24)
    set(h3, 'name',sprintf(strcat('cest_mtr','_',pulse_method,'_',sampling)),'numbertitle','off')
    
    
    Za(:,length(pools)+1)=Za;
    %each single pool except for free H2O with conc and exch rate set to 0
    for ii=2:length(pools)
        pools=pools_copy;
        pools{ii}{5}=0;
        pools{ii}{6}=0;
        [Za(:,ii), offset_w]=cest_sim_fast(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
            fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
        Za(:,ii)=Za(:,ii)/100;
        Za(:,ii)=Za(:,ii)/Za_norm*100;
    end
    
    %free H2O
    pools=pools_copy;
    for ii=2:length(pools)
        pools{ii}{5}=0;
        pools{ii}{6}=0;
    end
    [Za(:,1), offset_w]=cest_sim_fast(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
    Za(:,1)=Za(:,1)/100;
    Za(:,1)=Za(:,1)/Za_norm*100;
  
    
    
    %pure signal for each pool
    pureZa(:,1)=Za(:,1);
    
    for ii=2:length(pools)
        pureZa(:,ii)=100-(Za(:,ii)-Za(:,end));
    end
    
    figure(h3)
    subplot(1,2,1)
    for ii=1:length(pools)
        plot(offset_w,pureZa(:,ii),'Color',colors(ii,:),'LineWidth',4)
        hold on
    end
    legend(legend_pools)
end
end


function varargout=cest_sim_evolution(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
    fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par)

if strcmp(pulse_method,'ss')==1
    [Za, time_all, Za_all, sampling_rate_sat, sampling_rate_ex, offset, offset_order]=cest_sim_slave_ss_evolution(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
elseif strcmp(pulse_method,'ps')==1
    [Za, time_all, Za_all, sampling_rate_sat, sampling_rate_ex, offset, offset_order]=cest_sim_slave_ps_evolution(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
end

varargout{1}=Za;
varargout{2}=time_all;
varargout{3}=Za_all;
varargout{4}=sampling_rate_sat;
varargout{5}=sampling_rate_ex;
varargout{6}=offset;
varargout{7}=offset_order;
end


function [Za, offset_w]=cest_sim_fast(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
    fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par)

if strcmp(pulse_method,'ss')==1
    [Za, offset_w]=cest_sim_slave_ss_fast(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
elseif strcmp(pulse_method,'ps')==1
    [Za, offset_w]=cest_sim_slave_ps_fast(spacing,pools, offset_w, pulse_method, cwpe_approximation, sampling, ...
        fs, B1amp, pulse_dur, pulse_int, spoiler,pulse_repeat, dummy, recovery, pulse_sampling_factor,readout_par,other_par);
end
end



