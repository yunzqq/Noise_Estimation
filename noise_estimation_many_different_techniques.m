%%
% PhD Project, 2015-2019 
% Dionelis Nikolaos, CID: 00690438

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Clear all the previous data
clear all; clc; 
close all;

% Add the datapath of the voicebox
addpath ./voicebox
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

addpath ./clean

addpath ./babble/0dB
addpath ./babble/5dB
addpath ./babble/10dB
addpath ./babble/15dB

addpath ./airport/0dB
addpath ./airport/5dB
addpath ./airport/10dB
addpath ./airport/15dB

addpath ./car/0dB
addpath ./car/5dB
addpath ./car/10dB
addpath ./car/15dB
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapaths of the functions for Noise Estimation
addpath ./Functions_for_Noise_Estimation

% Add the datapaths of the functions for the Enhancement Systems
addpath ./Functions_for_Enhancement_System
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

addpath ./FCJF0
addpath ./NatoNoise0
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
filename = 'sp01_babble_sn15';
[s,fs] = readwav(fullfile(filename));

%s = resample(s,2,1);
s = v_resample(s,2,1);
fs = fs * 2;

z = s;
% z = [z; z];

% start_point1 = length(s);

filename = 'sp01';
[s2,fs2] = readwav(fullfile(filename));

s2 = v_resample(s2,2,1);
fs2 = fs2 * 2;

rrrqae = s - s2;
n = rrrqae;
% rrrqae = [rrrqae; rrrqae];

snr = 15;
z = v_addnoise(s2,fs,snr,'',n,fs);

% filename = 'sp01_babble_sn5';
% [s,fs] = readwav(fullfile(filename));
% 
% s = v_resample(s,2,1);
% fs = fs * 2;
% 
% z = [z; s; s];
% 
% filename = 'sp01';
% [s2,fs2] = readwav(fullfile(filename));
% 
% s2 = v_resample(s2,2,1);
% fs2 = fs2 * 2;
% 
% rrrqae2 = s - s2;
% rrrqae = [rrrqae; rrrqae2; rrrqae2];

figure; 
set(gcf,'Color','w');
subplot(2,2,1);
set(gca,'FontSize',15);
spgrambw(z,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noisy Speech Signal. 10 dB and 0 dB SNR. Babble NOIZEUS Noise.');

set(gcf,'Color','w');
subplot(2,2,2);
set(gca,'FontSize',15);
spgrambw(rrrqae,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noise Signal. Babble NOIZEUS Noise.');

set(gcf,'Color','w');
subplot(2,2,3);
set(gca,'FontSize',15);
spgrambw([s2;s2;s2],fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Clean Speech Signal. No Noise.');

true_ppx = true25_noise_parameters_my_proposed_algorithm( rrrqae );

[ppx, pvdfasf, signal_main_in_frames] = new_main2_main_proposed_algorithm_project( 10, z, rrrqae );
% pvdfasf is the SNR estimates

% lev = activlev(z,fs, 'd');
% %lev = lev * (10/38.6904);
% 
% lev2 = activlevg(z,fs, 'd');
% %lev2 = lev2 * (10/38.6904);

% pvdfasf is the SNR estimates
pvdfasf2 = mean(pvdfasf,2);

lev_total_in_every_frame = [];

sig_to_use = [];

for i = 1 : length(pvdfasf2)
   sig_to_use2 = signal_main_in_frames(i,:); 

   for i2 = 1 : length(pvdfasf2)
       sig_to_use = [sig_to_use sig_to_use2];
   end
   
   if (pvdfasf2(i) > 4) 
       % use P56
       lev = activlev(sig_to_use, fs, 'd');
       
   elseif (pvdfasf2(i) < -2)
       % use the harmonics summation principle
       lev = activlevg(sig_to_use, fs, 'd');
       
   else
       % use addition and smoothing
       % smoothing depends on the SNR
       if (-2 < pvdfasf2(i) <-1)
           ro_parameter = 0.08;
       elseif (-1 < pvdfasf2(i) <0)
           ro_parameter = 0.16;
       elseif (0 < pvdfasf2(i) <1)
           ro_parameter = 0.28;
       elseif (1 < pvdfasf2(i) <2)
           ro_parameter = 0.44;
       elseif (2 < pvdfasf2(i) <3)
           ro_parameter = 0.68;
       elseif (3 < pvdfasf2(i) <4)
           ro_parameter = 0.89;
       end
       
       lev2 = activlev(sig_to_use, fs, 'd');
       lev3 = activlevg(sig_to_use, fs, 'd');
       
       lev = (ro_parameter * lev2) + ((1-ro_parameter) * lev3);
       
   end
   
   lev_total_in_every_frame = [lev_total_in_every_frame lev];
end

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
ninc=0.010*fs;           % Frame increment for BW=200 Hz (in samples)
nwin=2*ninc;              % Frame length (in samples)
win=hamming(nwin);        % Analysis window
k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
%sf=abs(rfft(enframe(s,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array  

[afasfa sadfasfas] = size(ppx);
ppx55 = ppx(:,1:sadfasfas/2+1);

gsadgasdgsagdas = ppx55;
clear ppx55;

[afasfa sadfasfas] = size(gsadgasdgsagdas);
for loop62 = 1 : afasfa
    for loop52 = 1 : sadfasfas/2
        ppx55(loop62,loop52) = mean(gsadgasdgsagdas(loop62, (loop52-1)*2+1:(loop52*2)));
    end
end

subplot(2,2,4);
%figure;
set(gcf,'Color','w');
set(gca,'FontSize',15);
spgrambw(ppx55, [fs/ninc 0.5*(nwin+1)/fs fs/nwin], 'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Estimated Noise Signal. Babble NOIZEUS Noise.');

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% filename2 = 'SA2.WAV';
% [s,fs] = readsph((filename2));
% 
% % Read the noise file
% filename = 'factory1'; 
% [s2,fs2] = readwav((filename));
%
% % Define the SNR
% snr = 10;
%
% % Add noise from column vector n() at sample frequency fn with random start
% % Sample and wrapping around as required with a vorbis cross-fade
% z = v_addnoise(s,fs,snr,'',s2,fs2); 
% 
% zsafasfas = v_addnoise(s,fs,snr,'x',s2,fs2); 
% zsafasfas = z-zsafasfas(:,1);
% zsafasfas = zsafasfas(:,2);

%z = [z; z; z; z; z];
%rrrqae = [zsafasfas; zsafasfas; zsafasfas; zsafasfas; zsafasfas];
%rrrqae = zsafasfas;

filename = 'sp01_babble_sn10';
[s,fs] = readwav(fullfile(filename));

z = s;
% z = [z; z];

% start_point1 = length(s);

filename = 'sp01_16kHz';
[s2,fs2] = readwav(fullfile(filename));

rrrqae = s - s2;
% rrrqae = [rrrqae; rrrqae];

filename = 'sp01_babble_sn0';
[s,fs] = readwav(fullfile(filename));

z = [z; s; s];

filename = 'sp01_16kHz';
[s2,fs2] = readwav(fullfile(filename));

rrrqae2 = s - s2;
rrrqae = [rrrqae; rrrqae2; rrrqae2];

%xfinal = specsub_ns_ns_new(z, 'martin');

ninc=round(0.010*fs);   % frame increment [fs=sample frequency]
ovf=2;                  % overlap factor
f=rfft(enframe(z,hanning(ovf*ninc,'periodic'),ninc),2*ovf*ninc,2);
f=f.*conj(f);           % convert to power spectrum
xafafdaf = estnoiseg(f, ninc/fs);

[ggh ~] = size(xafafdaf);
for i = 1 : ggh
    prr(i,:) = [xafafdaf(i,:) flipud(xafafdaf(i,2:end-1))];
end

xafafdaf = prr;
clear prr;

xfinal = specsub_ns_ns_new2_new(z,xafafdaf);

% Store the clean speech signal
fs = 16000;
gg = [s2;s2;s2];
writewav(gg, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
result_fromSS = xfinal;
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final4 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
T = table(final4);
T.Properties.VariableNames = {'TestResults'};

% print the results
disp(['The following table shows the results that we just calculated.']);
disp(T);
disp(['   ']); disp(['   ']);

% set the figure data
f = figure; set(gcf,'Color','w');
set(f, 'Position', [500 400 300 100]);
dat =  {T{1,1};...
        T{2,1};...   
        T{3,1};...
        T{4,1};...
        T{5,1};...
        T{6,1};...
        T{7,1};...
        T{8,1};...
        T{9,1};...
        T{10,1};...
        T{11,1};...
        T{12,1};...
        T{13,1};...
        T{14,1};...
        T{15,1};...
        T{16,1};...
        T{17,1}; ...
        T{18,1}; ...
        T{19,1}; ...
        T{20,1}; ...
        T{21,1};};

% set the table as a figure
columnname =   {'Results'};
rnames = {'SNR mean', 'Segmental SNR', 'WSS mean', 'LLR mean', 'IS mean', 'CEP mean', ...
    'fw SNR segmental', 'SIG' 'BAK', 'OVL', 'SIG2' 'BAK2', 'OVL2', 'raw PESQ', 'Csig', 'Cbak', 'Covl', 'NCM' , 'CSh', 'CSm', 'CSl'};
columnformat = {'numeric'}; 
t = uitable('Units','normalized','Position',...
            [0.05 0.05 0.755 0.87], 'Data', dat,... 
            'ColumnName', columnname,...
            'RowName', rnames, ...
            'ColumnFormat', columnformat); 
        
%set the title
%set(gca,'FontSize',15); 
title('Table: MMSE Algorithm. The SNR falls in the signal. Only one speaker and one noise is tested.');        

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
filename = 'sp01_babble_sn10';
[s,fs] = readwav(fullfile(filename));

z = s;
% z = [z; z];

% start_point1 = length(s);

filename = 'sp01_16kHz';
[s2,fs2] = readwav(fullfile(filename));

rrrqae = s - s2;
% rrrqae = [rrrqae; rrrqae];

filename = 'sp01_babble_sn0';
[s,fs] = readwav(fullfile(filename));

z = [z; s; s];

filename = 'sp01_16kHz';
[s2,fs2] = readwav(fullfile(filename));

rrrqae2 = s - s2;
rrrqae = [rrrqae; rrrqae2; rrrqae2];

%xfinal = specsub_ns_ns_new(z, 'martin');

ninc=round(0.010*fs);   % frame increment [fs=sample frequency]
ovf=2;                  % overlap factor
f=rfft(enframe(z,hanning(ovf*ninc,'periodic'),ninc),2*ovf*ninc,2);
f=f.*conj(f);           % convert to power spectrum
xafafdaf = estnoiseg(f, ninc/fs);

[ggh ~] = size(xafafdaf);
for i = 1 : ggh
    prr(i,:) = [xafafdaf(i,:) flipud(xafafdaf(i,2:end-1))];
end

xafafdaf = prr;
clear prr;

xfinal = specsub_ns_ns_new2_new(z,xafafdaf);

figure; 
set(gcf,'Color','w');
subplot(2,2,1);
set(gca,'FontSize',15);
spgrambw(z,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noisy Speech Signal. 10 dB and 0 dB SNR. Babble NOIZEUS Noise.');

set(gcf,'Color','w');
subplot(2,2,2);
set(gca,'FontSize',15);
spgrambw(rrrqae,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noise Signal. Babble NOIZEUS Noise.');

set(gcf,'Color','w');
subplot(2,2,3);
set(gca,'FontSize',15);
spgrambw([s2;s2;s2],fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Clean Speech Signal. No Noise.');

subplot(2,2,4);
set(gcf,'Color','w');
set(gca,'FontSize',15);
spgrambw(xfinal, 16000, 'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Estimated Clean Speech Signal. MMSE Gerkmann. Babble NOIZEUS Noise.');

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%

% method = 'mcra2';
% [mcra2_noise_ps] = new_new2_noise_parameters(z, method);
% 
% true_ppx = true25_noise_parameters_my_proposed_algorithm( rrrqae );
% 
% method = 'conn_freq';
% [conn_freq_noise_ps] = new_new2_noise_parameters(z, method);
% 
% method = 'imcra';
% [imcra_noise_ps] = new_new2_noise_parameters(z, method);
% 
% method = 'mcra';
% [mcra_noise_ps] = new_new2_noise_parameters(z, method);
% 
% method = 'martin';
% [martin_noise_ps] = new_new2_noise_parameters(z, method);
% 
% method = 'doblinger';
% [doblinger_noise_ps] = new_new2_noise_parameters(z, method);
% 
% method = 'hirsch';
% [hirsch_noise_ps] = new_new2_noise_parameters(z, method);

filename = 'sp01_station_sn10';
[s,fs] = readwav(fullfile(filename));

z = s;
% z = [z; z];

% start_point1 = length(s);

filename = 'sp01_16kHz';
[s2,fs2] = readwav(fullfile(filename));

rrrqae = s - s2;
% rrrqae = [rrrqae; rrrqae];

filename = 'sp01_station_sn5';
[s,fs] = readwav(fullfile(filename));

z = [z; s; s];

filename = 'sp01_16kHz';
[s2,fs2] = readwav(fullfile(filename));

rrrqae2 = s - s2;
rrrqae = [rrrqae; rrrqae2; rrrqae2];

figure; 
set(gcf,'Color','w');
subplot(2,2,1);
set(gca,'FontSize',15);
spgrambw(z,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noisy Speech Signal. 10 dB and 0 dB SNR. Babble NOIZEUS Noise.');

set(gcf,'Color','w');
subplot(2,2,2);
set(gca,'FontSize',15);
spgrambw(rrrqae,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noise Signal. Babble NOIZEUS Noise.');

set(gcf,'Color','w');
subplot(2,2,3);
set(gca,'FontSize',15);
spgrambw([s2;s2;s2],fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Clean Speech Signal. No Noise.');

true_ppx = true25_noise_parameters_my_proposed_algorithm( rrrqae );

[ppx, pvdfasf] = new_main2_main_proposed_algorithm_project( 8, z, rrrqae );

ninc=0.010*fs;           % Frame increment for BW=200 Hz (in samples)
nwin=2*ninc;              % Frame length (in samples)
win=hamming(nwin);        % Analysis window
k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
%sf=abs(rfft(enframe(s,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array  

[afasfa sadfasfas] = size(ppx);
ppx55 = ppx(:,1:sadfasfas/2+1);

gsadgasdgsagdas = ppx55;
clear ppx55;

[afasfa sadfasfas] = size(gsadgasdgsagdas);
for loop62 = 1 : afasfa
    for loop52 = 1 : sadfasfas/2
        ppx55(loop62,loop52) = mean(gsadgasdgsagdas(loop62, (loop52-1)*2+1:(loop52*2)));
    end
end

subplot(2,2,4);
%figure;
set(gcf,'Color','w');
set(gca,'FontSize',15);
spgrambw(ppx55, [fs/ninc 0.5*(nwin+1)/fs fs/nwin], 'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Estimated Noise Signal. Station NOIZEUS Noise.');

% xfinal = SS_main_overSS_specsub_ns(z);

% ninc=0.010*fs;           % Frame increment for BW=200 Hz (in samples)
% nwin=2*ninc;              % Frame length (in samples)
% win=hamming(nwin);        % Analysis window
% k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
% %sf=abs(rfft(enframe(s,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array  
% 
% [afasfa sadfasfas] = size(ppx);
% ppx55 = ppx(:,1:sadfasfas/2+1);

% gsadgasdgsagdas = ppx55;
% clear ppx55;
% 
% [afasfa sadfasfas] = size(gsadgasdgsagdas);
% for loop62 = 1 : afasfa
%     for loop52 = 1 : sadfasfas/2
%         ppx55(loop62,loop52) = mean(gsadgasdgsagdas(loop62, (loop52-1)*2+1:(loop52*2)));
%     end
% end

subplot(2,2,4);
set(gcf,'Color','w');
set(gca,'FontSize',15);
spgrambw(xfinal,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Estimated Noise Signal. Babble NOIZEUS Noise.');

true_ppx = true25_noise_parameters_my_proposed_algorithm( rrrqae );

[ppx, pvdfasf] = new_main2_main_proposed_algorithm_project( 7, z, rrrqae );

[yp ~] = size(ppx);

%     fs = 16000;
%     ninc=round(0.020*fs/2);   % frame increment [fs=sample frequency]
%     ovf=2;                  % overlap factor
%     f=rfft(enframe(z,hanning(ovf*ninc,'periodic'),ninc),ovf*ninc,2);
%     f=f.*conj(f);           % convert to power spectrum
% 
%     %x=estnoiseg(f,ninc/fs); % estimate the noise power spectrum
%     % total_actbuf = estnoisem2(f,ninc/fs, ovf*ninc);
%     % tt = estnoisem3(total_actbuf, f,ninc/fs, ovf*ninc);
% 
%     % use transient1 
%     [total_actbuf,~,fffaa] = estnoisem(f, ninc/fs, ovf*ninc);
%     tt = estnoisem4(total_actbuf, fffaa, f, ninc/fs, ovf*ninc);
%     %tt = estnoisem(f, ninc/fs, ovf*ninc);
% 
%     [ao ~] = size(tt);
%     for i = 1 : ao
%         tt23(i,:) = [tt(i,:) flipud(tt(i,2:end-1))];
%     end
%     tt = tt23;
% 
%     future_tt = tt;
    
atasgas = [];
atasgas2 = [];
atasgas3 = [];
atasgas4 = [];
atasgas5 = [];
atasgas6 = [];
atasgas7 = [];
atasgas8 = [];
atasgas9 = [];
atasgas10 = [];

for i = 1 : yp
    if (i == 1) 
        var = 0.030;
    end
    
%     atasgas = [atasgas 10*log10(mcra2_noise_ps(i,20))]; 
%     atasgas2 = [atasgas2 10*log10(conn_freq_noise_ps(i,20))]; 
%     atasgas3 = [atasgas3 10*log10(imcra_noise_ps(i,20))]; 
%     atasgas4 = [atasgas4 10*log10(mcra_noise_ps(i,20))]; 
%     atasgas5 = [atasgas5 10*log10(martin_noise_ps(i,20))]; 
%     atasgas6 = [atasgas6 10*log10(doblinger_noise_ps(i,20))]; 
%     atasgas7 = [atasgas7 10*log10(hirsch_noise_ps(i,20))]; 

    atasgas8 = [atasgas8 10*log10(true_ppx(i,20))]; 
    
%     atasgas9 = [atasgas9 10*log10(future_tt(i,20))]; 
    
    atasgas10 = [atasgas10 10*log10(ppx(i,20))]; 

    var = var + (0.020/2);
end

% figure; set(gcf,'Color','w');
% %subplot(2,1,1);
% 
% set(gca,'FontSize',15);
% plot(0.020:(0.030/2):var-(0.030/2), atasgas, '--*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas2, '--r*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas3, '--k*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas4, '--c*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas5, '--m*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas6, ':o');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas7, ':ro');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas8, ':ok');
% hold on;
% 
% %plot(0.020:(0.030/2):var-(0.030/2), atasgas9, 'o');
% %hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas10, ':go');
% hold off;
% 
% set(gca,'FontSize',15);
% title('Babble Noise.  SNR falls at t = 4.2 s.   Noise Estimation Algorithms.   20 ms with 10 ms overlap. Power (in dB) vs Time (s) for one frequency bin.');
% 
% xlabel('Time in seconds (s)');
% ylabel('Power in dB');
% 
% legend(' MCRA2',' Connected TF', ' IMCRA', ' MCRA', ' Martin MS', ' Doblinger', ' Hirsch', ' True', ' Proposed, 1.1', 'Location', 'SouthEast');
% 
% grid on; grid minor;

% ninc=0.010*fs;           % Frame increment for BW=200 Hz (in samples)
% nwin=2*ninc;              % Frame length (in samples)
% win=hamming(nwin);        % Analysis window
% k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
% %sf=abs(rfft(enframe(s,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array  
% 
% [afasfa sadfasfas] = size(ppx);
% ppx55 = ppx(:,1:sadfasfas/2+1);
% 
% % gsadgasdgsagdas = ppx55;
% % clear ppx55;
% % 
% % [afasfa sadfasfas] = size(gsadgasdgsagdas);
% % for loop62 = 1 : afasfa
% %     for loop52 = 1 : sadfasfas/2
% %         ppx55(loop62,loop52) = mean(gsadgasdgsagdas(loop62, (loop52-1)*2+1:(loop52*2)));
% %     end
% % end
% 
% subplot(2,2,4);
% set(gcf,'Color','w');
% set(gca,'FontSize',15);
% spgrambw(ppx55,[fs/ninc 0.5*(nwin+1)/fs fs/nwin],'pJcw');
% set(gca,'FontSize',15);
% title('Spectrogram. Estimated Noise Signal. Leopard Nato Noise.');

figure; set(gcf,'Color','w');
subplot(2,1,1);
set(gca,'FontSize',15);
plot(0.020:(0.020/2):var-(0.020/2), atasgas8, 'r');
hold on;

plot(0.020:(0.020/2):var-(0.020/2), real(atasgas10));

% hold on;
% plot(0.020:(0.020/2):var-(0.020/2), 10*log10(true_ppx2(:,20)), 'g');
hold off;

set(gca,'FontSize',15);
title('Babble Noise.  SNR falls at t = 4.2 s.   Estimated Noise.  20 ms with 10 ms overlap. Noise Power (in dB) vs Time (s) for one frequency bin.');

xlabel('Time in seconds (s)');
ylabel('Power in dB');

legend(' True, No Smoothing', ' Proposed, 2.1, non-voiced when pv < 0.2', ' True, Smoothing', 'Location', 'SouthEast');

grid on; grid minor;

subplot(2,1,2);
plot(0.020:(0.020/2):var-(0.020/2), pvdfasf);
set(gca,'FontSize',15);
title('Frame Voiced Probability pv.');

xlabel('Time in seconds (s)');
ylabel('pv');
grid on; grid minor;


% figure;
% plot(0.020:(0.030/2):var-(0.030/2), 10.^(atasgas8/10), 'r');
% hold on;
% plot(0.020:(0.030/2):var-(0.030/2), 10.^(real(atasgas10)/10));
% hold off;
% grid on; grid minor;

% 
% atasgas = [];
% atasgas2 = [];
% atasgas3 = [];
% atasgas4 = [];
% atasgas5 = [];
% atasgas6 = [];
% atasgas7 = [];
% atasgas8 = [];
% atasgas9 = [];
% atasgas10 = [];
% 
% for i = 1 : yp
%     if (i == 1) 
%         var = 0.030;
%     end
%     
%     atasgas = [atasgas 10*log10(mcra2_noise_ps(i,300))]; 
%     atasgas2 = [atasgas2 10*log10(conn_freq_noise_ps(i,300))]; 
%     atasgas3 = [atasgas3 10*log10(imcra_noise_ps(i,300))]; 
%     atasgas4 = [atasgas4 10*log10(mcra_noise_ps(i,300))]; 
%     atasgas5 = [atasgas5 10*log10(martin_noise_ps(i,300))]; 
%     atasgas6 = [atasgas6 10*log10(doblinger_noise_ps(i,300))]; 
%     atasgas7 = [atasgas7 10*log10(hirsch_noise_ps(i,300))]; 
%     atasgas8 = [atasgas8 10*log10(true_ppx(i,300))]; 
%     atasgas9 = [atasgas9 10*log10(future_tt(i,300))]; 
%     atasgas10 = [atasgas10 10*log10(ppx(i,300))]; 
% 
%     var = var + (0.030/2);
% end
% 
% figure; set(gcf,'Color','w');
% %subplot(2,1,1);
% 
% set(gca,'FontSize',15);
% plot(0.020:(0.030/2):var-(0.030/2), atasgas, '--*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas2, '--r*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas3, '--k*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas4, '--c*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas5, '--m*');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas6, ':o');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas7, ':ro');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas8, ':ok');
% hold on;
% 
% %plot(0.020:(0.030/2):var-(0.030/2), atasgas9, 'o');
% %hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), atasgas10, ':go');
% hold off;
% 
% set(gca,'FontSize',15);
% title('Babble Noise.  SNR falls at t = 4.2 s.   Noise Estimation Algorithms.   20 ms with 10 ms overlap. Power (in dB) vs Time (s) for one frequency bin.');
% 
% xlabel('Time in seconds (s)');
% ylabel('Power in dB');
% 
% legend(' MCRA2',' Connected TF', ' IMCRA', ' MCRA', ' Martin MS', ' Doblinger', ' Hirsch', ' True', ' Proposed, 1.1', 'Location', 'SouthEast');
% 
% grid on; grid minor;
% 
% figure; set(gcf,'Color','w');
% subplot(2,1,1);
% set(gca,'FontSize',15);
% plot(0.020:(0.030/2):var-(0.030/2), atasgas8, 'r');
% hold on;
% 
% plot(0.020:(0.030/2):var-(0.030/2), real(atasgas10));
% hold off;
% 
% set(gca,'FontSize',15);
% title('Babble Noise.  SNR falls at t = 4.2 s.   Noise Estimation Algorithms.   20 ms with 10 ms overlap. Power (in dB) vs Time (s) for one frequency bin.');
% 
% xlabel('Time in seconds (s)');
% ylabel('Power in dB');
% 
% legend(' True', ' Proposed, 1.1', 'Location', 'SouthEast');
% 
% grid on; grid minor;
% 
% subplot(2,1,2);
% plot(0.020:(0.030/2):var-(0.030/2), pvdfasf);
% set(gca,'FontSize',15);
% title('Voiced Frame Probability pv.');
% 
% xlabel('Time in seconds (s)');
% ylabel('pv');
% grid on; grid minor;
% 
% 
% figure;
% plot(0.020:(0.030/2):var-(0.030/2), 10.^(atasgas8/10), 'r');
% hold on;
% plot(0.020:(0.030/2):var-(0.030/2), 10.^(real(atasgas10)/10));
% hold off;
% grid on; grid minor;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%

% filename = 'sp01_babble_sn5';
% [s,fs] = readwav(fullfile(filename));

    fs = 16000;
    ninc=round(0.020*fs/2);   % frame increment [fs=sample frequency]
    ovf=2;                  % overlap factor
    f=rfft(enframe(z,hanning(ovf*ninc,'periodic'),ninc),ovf*ninc,2);
    f=f.*conj(f);           % convert to power spectrum

    %x=estnoiseg(f,ninc/fs); % estimate the noise power spectrum
    % total_actbuf = estnoisem2(f,ninc/fs, ovf*ninc);
    % tt = estnoisem3(total_actbuf, f,ninc/fs, ovf*ninc);

    % use transient1 
    [total_actbuf,~,fffaa] = estnoisem(f, ninc/fs, ovf*ninc);
    tt = estnoisem4(total_actbuf, fffaa, f, ninc/fs, ovf*ninc);
    %tt = estnoisem(f, ninc/fs, ovf*ninc);

    [ao ~] = size(tt);
    for i = 1 : ao
        tt23(i,:) = [tt(i,:) flipud(tt(i,2:end-1))];
    end
    tt = tt23;

    future_tt = tt;
    ppx = future_tt;

%ppx = new_main2_main_proposed_algorithm_project( 9, z, rrrqae );

% with smoothing
true_ppx = true2_noise_parameters_my_proposed_algorithm( rrrqae );

[MSE, RMS_error] = MSE_function(true_ppx, ppx);

[~, MedSE_final] = MedSE_function(true_ppx, ppx);

[std_final, iqr_final] = Std_function(ppx);

tt = ppx;
true_ppp_total_noise_ps = true_ppx;

figure; 
set(gcf,'Color','w');
subplot(2,2,1);
set(gca,'FontSize',15);
spgrambw(z,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noisy Speech Signal. 5 dB SNR. Leopard Nato Noise.');

set(gcf,'Color','w');
subplot(2,2,2);
set(gca,'FontSize',15);
spgrambw(rrrqae,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noise Signal. Leopard Nato Noise.');

set(gcf,'Color','w');
subplot(2,2,3);
set(gca,'FontSize',15);
spgrambw(s,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Clean Speech Signal. No Noise.');

ninc=0.010*fs;           % Frame increment for BW=200 Hz (in samples)
nwin=2*ninc;              % Frame length (in samples)
win=hamming(nwin);        % Analysis window
k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
%sf=abs(rfft(enframe(s,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array  

[afasfa sadfasfas] = size(ppx);
ppx55 = ppx(:,1:sadfasfas/2+1);

% gsadgasdgsagdas = ppx55;
% clear ppx55;
% 
% [afasfa sadfasfas] = size(gsadgasdgsagdas);
% for loop62 = 1 : afasfa
%     for loop52 = 1 : sadfasfas/2
%         ppx55(loop62,loop52) = mean(gsadgasdgsagdas(loop62, (loop52-1)*2+1:(loop52*2)));
%     end
% end

subplot(2,2,4);
set(gcf,'Color','w');
set(gca,'FontSize',15);
spgrambw(ppx55,[fs/ninc 0.5*(nwin+1)/fs fs/nwin],'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Estimated Noise Signal. Leopard Nato Noise.');

% graph 2
% power vs time
[yp ~] = size(tt);

atasgas = [];
atasgas2 = [];

for i = 1 : yp
    if (i == 1) 
        var = 0.020;
    end
    
    atasgas = [atasgas 10*log10(tt(i,20))]; 
    atasgas2 = [atasgas2 10*log10(true_ppp_total_noise_ps(i,20))]; 

    var = var + (0.020/2);
end

figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
plot(0.020:(0.020/2):var, atasgas);
%plot(0.020:(0.020/2):var-(0.020/2), atasgas2+mean(atasgas2)*randn(size(atasgas2)), 'r');
hold on;

%plot(0.020:(0.020/2):var-(0.020/2), atasgas2);
plot(0.020:(0.020/2):var, atasgas2, 'r');
hold off;

set(gca,'FontSize',15);
title('0 dB SNR. 20 ms with 10 ms overlap.  Leopard Noise. Nato Database.  Power (in dB) vs Time for one frequency bin. This is for one frequency bin.');
xlabel('Time in seconds (s)');
ylabel('Power in dB');

legend(' Estimated Power, Proposed 2.2',' True Power, With Smoothing');
grid on; grid minor;

% filename = 'sp01_babble_sn5';
% [s,fs] = readwav(fullfile(filename));
% 
% z = s;
% z = [z; z];
% 
% filename = 'sp01_16kHz';
% [s2,fs2] = readwav(fullfile(filename));
% 
% rrrqae = s - s2;
% % rrrqae = zsafasfas;
% 
% rrrqae = [rrrqae; rrrqae];
% rrrqae2 = rrrqae;
% 
% filename = 'sp01_babble_sn5';
% [s,fs] = readwav(fullfile(filename));
% 
% z2 = z;
% z = s;
% z = [z2; z; z; z];
% 
% filename = 'sp01_16kHz';
% [s2,fs2] = readwav(fullfile(filename));
% 
% rrrqae = s - s2;
% % rrrqae = zsafasfas;
% 
% rrrqae = [rrrqae2; rrrqae; rrrqae; rrrqae];

fs = 16000;

% figure; 
% set(gcf,'Color','w');
% subplot(2,2,1);
% set(gca,'FontSize',15);
% spgrambw(z,fs,'pJcw');
% set(gca,'FontSize',15);
% title('Spectrogram. Noisy Speech Signal. 5 dB SNR. Station Noise.');
% 
% set(gcf,'Color','w');
% subplot(2,2,2);
% set(gca,'FontSize',15);
% spgrambw(rrrqae,fs,'pJcw');
% set(gca,'FontSize',15);
% title('Spectrogram. Noise Signal. Station Noise.');
% 
% set(gcf,'Color','w');
% subplot(2,2,3);
% set(gca,'FontSize',15);
% spgrambw([s2;s2;s2;s2;s2],fs,'pJcw');
% set(gca,'FontSize',15);
% title('Spectrogram. Clean Speech Signal. No Noise.');

%writewav(z, 16000, '1asfdasfdasaasfa');
%tt = new_noise_parameters_my_proposed_algorithm( '1asfdasfdasaasfa' );
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%

fs = 16000;
ninc=round(0.020*fs/2);   % frame increment [fs=sample frequency]
ovf=2;                  % overlap factor
f=rfft(enframe(z,hanning(ovf*ninc,'periodic'),ninc),2*ovf*ninc,2);
f=f.*conj(f);           % convert to power spectrum

%x=estnoiseg(f,ninc/fs); % estimate the noise power spectrum
% total_actbuf = estnoisem2(f,ninc/fs, ovf*ninc);
% tt = estnoisem3(total_actbuf, f,ninc/fs, ovf*ninc);

% use transient1 
[total_actbuf,~,fffaa] = estnoisem(f, ninc/fs, ovf*ninc);
tt = estnoisem4(total_actbuf, fffaa, f, ninc/fs, ovf*ninc);

ttasfttasdfas = estnoisem(f,ninc/fs, ovf*ninc);

[ao ~] = size(tt);
for i = 1 : ao
    tt23(i,:) = [tt(i,:) flipud(tt(i,2:end-1))];
end
tt = tt23;

future_tt = tt;

true_ppx = true2_noise_parameters_my_proposed_algorithm( rrrqae );

true_ppp_total_noise_ps = true_ppx;

% graph 2
% power vs time
[yp ~] = size(tt);

atasgas = [];
atasgas2 = [];
atasgas3 = [];

for i = 1 : yp
    if (i == 1) 
        var = 0.020;
    end
    
    atasgas = [atasgas 10*log10(tt(i,20))]; 
    atasgas2 = [atasgas2 10*log10(true_ppp_total_noise_ps(i,20))]; 
    atasgas3 = [atasgas3 10*log10(ttasfttasdfas(i,20))]; 
    
    var = var + (0.020/2);
end

figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
plot(0.020:(0.020/2):var, atasgas);
%plot(0.020:(0.020/2):var-(0.020/2), atasgas2+mean(atasgas2)*randn(size(atasgas2)), 'r');
hold on;

%plot(0.020:(0.020/2):var-(0.020/2), atasgas2);
plot(0.020:(0.020/2):var, atasgas2, 'r');

hold on;

%plot(0.020:(0.020/2):var-(0.020/2), atasgas2);
plot(0.020:(0.020/2):var, atasgas3, 'k');
hold off;

set(gca,'FontSize',15);
title('20 ms with 10 ms overlap.  SNR changes. Tracking Delay.  Babble Noise. Power (in dB) vs Time for one frequency bin.');
xlabel('Time in seconds (s)');
ylabel('Power in dB');

legend(' Estimated Power, Transient algorithm',' True Power, With Smoothing', ' Estimated Power, MS algorithm', 'Location', 'South');
grid on; grid minor;

ppx = new_main2_main_proposed_algorithm_project( 8, z, rrrqae );
%[ppx, true_ppx] = main_proposed_algorithm_project( 6, z, rrrqae );

% ppx = future_tt;

% with smoothing
true_ppx = true2_noise_parameters_my_proposed_algorithm( rrrqae );

[MSE, RMS_error] = MSE_function(true_ppx, ppx);

[~, MedSE_final] = MedSE_function(true_ppx, ppx);

[std_final, iqr_final] = Std_function(ppx);

% Frequency increment (in hertz per sample):
[ada pooii] = size(true_ppx);
N = pooii/2;
Fs = 16000;
dF = Fs/N;

% Frequency domain axes
% in hertz
fr = -(Fs/2):dF/2:Fs/2-dF/2; 

% graph 1
% power vs frequency for one frame
figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
plot(fr, 10*log10(ppx(100,:)));
hold on;
plot(fr, 10*log10(true_ppx(100,:)), 'r');
hold off;

set(gca,'FontSize',15);
title('Power (in dB) vs Frequency for one frame. This is one frame.  Buccaneer2 Noise. Nato Noise.  SNR = 5 dB.');
xlabel('Frequency in Hz');
ylabel('Power in dB');

grid on; grid minor;
legend(' Estimated Power, Proposed 2.2',' True Power, Smoothing was used');

tt = ppx;
true_ppp_total_noise_ps = true_ppx;

% graph 2
% power vs time
[yp ~] = size(tt);

atasgas = [];
atasgas2 = [];

for i = 1 : yp
    if (i == 1) 
        var = 0.020;
    end
    
    atasgas = [atasgas 10*log10(tt(i,20))]; 
    atasgas2 = [atasgas2 10*log10(true_ppp_total_noise_ps(i,20))]; 

    var = var + (0.020/2);
end

figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
plot(0.020:(0.020/2):var, atasgas);
%plot(0.020:(0.020/2):var-(0.020/2), atasgas2+mean(atasgas2)*randn(size(atasgas2)), 'r');
hold on;

%plot(0.020:(0.020/2):var-(0.020/2), atasgas2);
plot(0.020:(0.020/2):var, atasgas2, 'r');
hold off;

set(gca,'FontSize',15);
title('20 ms with 10 ms overlap. Volvo Car Noise. Power (in dB) vs Time for one frequency bin. This is for one frequency bin.');
xlabel('Time in seconds (s)');
ylabel('Power in dB');

legend(' Estimated Power, Proposed 2.2',' True Power, With Smoothing');
grid on; grid minor;
%, 'Location','southwest'

% ninc=0.010*fs;           % Frame increment for BW=200 Hz (in samples)
% nwin=2*ninc;              % Frame length (in samples)
% win=hamming(nwin);        % Analysis window
% k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
% %sf=abs(rfft(enframe(s,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array  
% 
% [afasfa sadfasfas ] = size(ppx);
% ppx55 = ppx(:,1:sadfasfas/2+1);
% 
% set(gcf,'Color','w');
% subplot(2,2,4);
% set(gca,'FontSize',15);
% spgrambw(ppx55,[fs/ninc 0.5*(nwin+1)/fs fs/nwin],'pJcw');
% set(gca,'FontSize',15);
% title('Spectrogram. Estimated Noise Signal. Station Noise.');

% [MSE, RMS_error] = MSE_function(true_ppx, ppx);

% Frequency increment (in hertz per sample):
[ada pooii] = size(true_ppx);
N = pooii/2;
Fs = 16000;
dF = Fs/N;

% Frequency domain axes
% in hertz
fr = -(Fs/2):dF/2:Fs/2-dF/2; 

% graph 1
% power vs frequency for one frame
figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
plot(fr, 10*log10(ppx(100,:)));
hold on;
plot(fr, 10*log10(true_ppx(100,:)), 'r');
hold off;

set(gca,'FontSize',15);
title('Power (in dB) vs Frequency for one frame. This is one frame.');
xlabel('Frequency in Hz');
ylabel('Power in dB');

grid on; grid minor;
legend(' Estimated Power, Proposed 2.1',' True Power, Smoothing was used');

tt = ppx;
true_ppp_total_noise_ps = true_ppx;

% graph 2
% power vs time
[yp ~] = size(tt);

atasgas = [];
atasgas2 = [];

for i = 1 : yp
    if (i == 1) 
        var = 0.020;
    end
    
    atasgas = [atasgas 10*log10(tt(i,100))]; 
    atasgas2 = [atasgas2 10*log10(true_ppp_total_noise_ps(i,100))]; 

    var = var + (0.020/2);
end

figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
plot(0.020:(0.020/2):var, atasgas);
%plot(0.020:(0.020/2):var-(0.020/2), atasgas2+mean(atasgas2)*randn(size(atasgas2)), 'r');
hold on;

%plot(0.020:(0.020/2):var-(0.020/2), atasgas2);
plot(0.020:(0.020/2):var, atasgas2, 'r');
hold off;

set(gca,'FontSize',15);
title('20 ms with 10 ms overlap.  SNR = 5 dB. Babble Noise. Power (in dB) vs Time for one frequency bin. This is for one frequency bin.');
xlabel('Time in seconds (s)');
ylabel('Power in dB');

legend(' Estimated Power, Proposed 1.1',' True Power, With Smoothing');
grid on; grid minor;

H = ppx;
[poiii pooii] = size(H);

w = 1:poiii;
% freq_axis = 0:pooii-1;

% Frequency increment (in hertz per sample):
N = pooii;
Fs = 16000;
dF = Fs/N;

% Frequency domain axes
% in hertz
fr = 0:dF/2:Fs/2-dF/2; 

% display the result
figure; 
subplot(2,2,3);
set(gcf,'Color','w');

% we use: mesh(X,Y,Z);
mesh(w, fr, H');

% set the title of the image
set(gca,'FontSize',15);
title('Results of the 3D plot. Estimated Power.');
zlabel('PSD estimate (dB)'); 
ylabel('Frequency Axis in Hz'); 
xlabel('Time index in frames (frame index)');

H = true_ppx;
[poiii pooii] = size(H);

w = 1:poiii;
% freq_axis = 0:pooii-1;

% Frequency increment (in hertz per sample):
N = pooii;
Fs = 16000;
dF = Fs/N;

% Frequency domain axes
% in hertz
fr = 0:dF/2:Fs/2-dF/2; 

% display the result
subplot(2,1,1);
set(gcf,'Color','w');

% we use: mesh(X,Y,Z);
mesh(w, fr, H');

% set the title of the image
set(gca,'FontSize',15);
title('Results of the 3D plot. True Power.');
zlabel('PSD estimate (dB)'); 
ylabel('Frequency Axis in Hz'); 
xlabel('Time index in frames (frame index)');

H = ppx_8;
[poiii pooii] = size(H);

w = 1:poiii;
% freq_axis = 0:pooii-1;

% Frequency increment (in hertz per sample):
N = pooii;
Fs = 16000;
dF = Fs/N;

% Frequency domain axes
% in hertz
fr = 0:dF/2:Fs/2-dF/2; 

% display the result
subplot(2,2,4);
set(gcf,'Color','w');

% we use: mesh(X,Y,Z);
mesh(w, fr, H');

% set the title of the image
set(gca,'FontSize',15);
title('Results of the 3D plot. Estimated Power.');
zlabel('PSD estimate (dB)'); 
ylabel('Frequency Axis in Hz'); 
xlabel('Time index in frames (frame index)');

ninc=(16000*0.020)/2;           

% Frame length (in samples)
nwin=16000*0.020;              

win=hamming(nwin);        % Analysis window
k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
sf=abs(rfft(enframe(input_main_noisy_speech_signal_z,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array                


spgrambw(sf,[fs/ninc 0.5*(nwin+1)/fs fs/nwin]);  % Plot spectrum array
spgrambw(sf,[fs/ninc 0.5*(nwin+1)/fs fs/nwin],'Jc');  % Plot spectrum array
ninc=(16000*0.020)/2;           

% Frame length (in samples)
nwin=16000*0.020;              

win=hamming(nwin);        % Analysis window
k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
sf=abs(rfft(enframe(x_final,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array                
spgrambw(sf,[fs/ninc 0.5*(nwin+1)/fs fs/nwin],'Jc');  % Plot spectrum array
Undefined function or variable 'x_final'.
 
Did you mean:
sf=abs(rfft(enframe(xfinal,win,ninc),nwin,2)).^2/k           ;
spgrambw(sf,[fs/ninc 0.5*(nwin+1)/fs fs/nwin],'Jc');  % Plot spectrum array
figure; set(gcf,'Color','w');
set(gca,'FontSize',15); spgrambw(sf,[fs/ninc 0.5*(nwin+1)/fs fs/nwin], 'pJcw');  % Plot spectrum array
title('Enhanced Speech Signal. After Processing. Station Noise. 5 dB SNR. ');

ninc=(16000*0.020)/2;           

% Frame length (in samples)
nwin=16000*0.020;              

win=hamming(nwin);        % Analysis window
k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
sf=abs(rfft(enframe(input_main_noisy_speech_signal_z,win,ninc),nwin,2)).^2/k; 
figure; set(gcf,'Color','w');
set(gca,'FontSize',15); spgrambw(sf,[fs/ninc 0.5*(nwin+1)/fs fs/nwin], 'pJcw');  % Plot spectrum array
title('Noisy Speech Signal. Before Processing. Station Noise. 5 dB SNR. ');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Frame increment for BW=200 Hz (in samples)
ninc=(16000*0.020)/2;           

% Frame length (in samples)
nwin=16000*0.020;              

win=hamming(nwin);        % Analysis window
k=0.5*fs*sum(win.^2);     % Scale factor to convert to power/Hz
sf=abs(rfft(enframe(z,win,ninc),nwin,2)).^2/k;           % Calculate spectrum array                


% spgrambw(sf,[fs/ninc 0.5*(nwin+1)/fs fs/nwin]);  % Plot spectrum array
% [fx,tx,pv,fv] = fxpefac(sf, 16000, [fs/ninc 0.5*(nwin+1)/fs fs/nwin]);
[fx,tx,pv,fv] = fxpefac(sf, [fs/ninc 0.5*(nwin+1)/fs fs/nwin]);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% snr = 5;
% 
% % Add noise from column vector n() at sample frequency fn with random start
% % Sample and wrapping around as required with a vorbis cross-fade
% z = v_addnoise(s,fs,snr,'',s2,fs2); 
% 
% z = [z; z; z; z; z];
% 
% zfasfa3 = v_addnoise(s,fs,snr,'x',s2,fs2); 
% rrrqae = zfasfa3(:,2);
% 
% rrrqae = [rrrqae; rrrqae; rrrqae; rrrqae; rrrqae];

% writewav(z, 16000, '1asfdasfdasaasfa');
% writewav(rrrqae, 16000, 'noise_only_1asfdasfdasaasfa');
% 
% fs = 16000;
% ninc=round(0.020*fs);   % frame increment [fs=sample frequency]
% ovf=2;                  % overlap factor
% f=rfft(enframe(z,hanning(ovf*ninc,'periodic'),ninc),ovf*ninc,2);
% f=f.*conj(f);           % convert to power spectrum
% 
% %x=estnoiseg(f,ninc/fs); % estimate the noise power spectrum
% tt = estnoisem(f,ninc/fs);
% 
% % true_ppp_total_noise_ps = true_noise_parameters_my_proposed_algorithm( rrrqae );
% 
% tt = new_noise_parameters_my_proposed_algorithm( '1asfdasfdasaasfa' );

[ tt, true_ppp_total_noise_ps ] = main_proposed_algorithm_project( 8, z, rrrqae );

%true_ppp_total_noise_ps = true2_noise_parameters_my_proposed_algorithm( 'noise_only_1asfdasfdasaasfa' );

% graph 1
% power vs frequency for one frame
figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
plot(10*log10(tt(100,:)));
hold on;
plot(10*log10(true_ppp_total_noise_ps(100,:)), 'r');
hold off;

set(gca,'FontSize',15);
title('Power (in dB) vs Frequency for one frame. This is one frame.');
xlabel('Frequency in Frequency Bins');
ylabel('Power in dB');

legend(' Estimated Power, Proposed',' True Power');

% graph 2
% power vs time

[yp ~] = size(tt);

atasgas = [];
atasgas2 = [];

for i = 1 : yp
    if (i == 1) 
        var = 0.020;
    end
    
    atasgas = [atasgas 10*log10(tt(i,81))]; 
    atasgas2 = [atasgas2 10*log10(true_ppp_total_noise_ps(i,81))]; 
    
    var = var + (0.020/2);
end

figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
plot(0.020:(0.020/2):var, atasgas);
%plot(0.020:(0.020/2):var-(0.020/2), atasgas2+mean(atasgas2)*randn(size(atasgas2)), 'r');
hold on;

%plot(0.020:(0.020/2):var-(0.020/2), atasgas2);
plot(0.020:(0.020/2):var, atasgas2, 'r');
hold off;

set(gca,'FontSize',15);
title('SNR = 5 dB. Power (in dB) vs Time for one frequency bin. This is for one frequency bin.');
xlabel('Time in seconds (s)');
ylabel('Power in dB');
legend(' Estimated Power, Proposed', ' True Power');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% find MSE

% We now use the "NOIZEUS_Noisy_Speech_Signal_Database_16kHz" database.
%filename = 'sp01_airport_sn0';
%filename = 'sp01_airport_sn5';

filename = 'sp01_babble_sn5';
%filename = 'sp01_airport_sn10';
%filename = 'sp01_airport_sn15';

MSE_3_total = [];
MSE_notRel_3_total = [];
MSE_med_3_total = [];
RMS1 = [];

% Outter for loop
%for outter_loop = 1 : 4
for outter_loop = 1 : 4
    %for inner_loop = 1 : 30
    for inner_loop = 1 : 30
        
%         if (inner_loop == 25)
%            inner_loop = 26; 
%         end
        
        if (outter_loop == 3 && inner_loop == 4)
            inner_loop = 5;
        end


        % We change the filename based on inner_loop
        % We use a temporary variable
        temp_variable = num2str(inner_loop);

        % We now change the filename based on inner_loop
        if (inner_loop < 10) 
            filename(3) = num2str(0);
            filename(4) = num2str(temp_variable(1));
        else
            filename(3) = num2str(temp_variable(1));
            filename(4) = num2str(temp_variable(2));
        end

        % If statement to change the filename
        if (outter_loop == 1) 
            %filename(16) = num2str(0);
            filename(15) = num2str(0);
        elseif (outter_loop == 2) 
            %filename(16) = num2str(5);
            filename(15) = num2str(5);
        elseif (outter_loop == 3) 
            %filename(16) = num2str(1);
            filename(15) = num2str(1);
            
            %filename(17) = num2str(0);
            filename(16) = num2str(0);
        else 
            %filename(16) = num2str(1);
            filename(15) = num2str(1);
            
            %filename(17) = num2str(5);
            filename(16) = num2str(5);
        end

        % Read the noisy speech file
        % [input_main_noisy_speech_signal_z, ~] = readwav((filename));
        
        filename99 = filename(1:4);
        
        % use: sp01_16kHz
        filename99(5) = '_';
        filename99(6) = num2str(1);
        filename99(7) = num2str(6);
        filename99(8) = 'k';
        filename99(9) = 'H';
        filename99(10) = 'z';

        % Read the noise file
        % [true, ~] = readwav((filename99));
        
        % input_noise_noiseOnly_signal = input_main_noisy_speech_signal_z - true;
        
        %[MSE_12, MSE_notRelative_12, MedSE_final_12] = main_proposed_algorithm_project( 8, input_main_noisy_speech_signal_z, input_noise_noiseOnly_signal );

        filename
        
        total_noise_ps = new_noise_parameters_my_proposed_algorithm( filename );

        true_total_noise_ps = true_noise_parameters_my_proposed_algorithm( filename, filename99 );

        [MSE_3, RMS_error_3] = MSE_function(true_total_noise_ps, total_noise_ps);
        [MSE_notRelative_3, MedSE_final_3] = MedSE_function(true_total_noise_ps, total_noise_ps);
        
        MSE_3_total = [MSE_3_total MSE_3];
        RMS1 = [RMS1 RMS_error_3];
        MSE_notRel_3_total = [MSE_notRel_3_total MSE_notRelative_3];
        MSE_med_3_total = [MSE_med_3_total MedSE_final_3];       
        
    end
end

MSE_3_final_total = mean(MSE_3_total);
final_MSE_notRel_3_total = mean(MSE_notRel_3_total);
final_MSE_med_3_total = mean(MSE_med_3_total);
RMS1_total = mean(RMS1); 
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% % Read the PHN TIMIT file
% filename2 = 'SA2.WAV';
% [s,fs] = readsph((filename2));
% 
% % Read the noise file
% filename = 'white'; 
% [s2,fs2] = readwav((filename));
% 
% % Define the SNR
% snr = 0;
% 
% % Add noise from column vector n() at sample frequency fn with random start
% % Sample and wrapping around as required with a vorbis cross-fade
% z = v_addnoise(s,fs,snr,'x',s2,fs2); 

filename = 'sp01_train_sn5';
s = readwav(fullfile(filename));
fs = 16000;

filename = 'sp01';
[s2,fs2] = readwav(fullfile(filename));
f = resample(s2,2,1);

% Define the noise signal
noise_signal = s-f;

% Define s
s = noise_signal;

fs = 16000;
figure; set(gcf,'Color','w');
set(gca,'FontSize',15);
spgrambw(s,fs,'pJcw');
set(gca,'FontSize',15);
title('Spectrogram. Noise only: train.');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% Plot noisy signals in the time domain

% Open the file
filename = 'sp01_train_sn0';
s = wavread(fullfile(filename));

% Plot the result
figure; set(gcf,'Color','w');
subplot(2,1,1);
plot(0:((45058/16000)/length(s)):(45058-((45058/length(s))))/16000,s)
set(gca,'FontSize',15);
title('SNR = 0 dB. Time domain representation of noisy speech signal.'); 
xlabel('Time in Seconds (seconds)');
ylabel('Amplitude');
grid on; grid minor;

% Open the file
filename = 'sp01_train_sn15';
s = wavread(fullfile(filename));

% Plot the result
subplot(2,1,2);
plot(0:((45058/16000)/length(s)):(45058-((45058/length(s))))/16000,s)
set(gca,'FontSize',15);
title('SNR = 15 dB. Time domain representation of noisy speech signal.');
xlabel('Time in Seconds (seconds)');
ylabel('Amplitude');
grid on; grid minor;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% We create the noisy signals.

% % Define the matrix for SNR 15 
% for i = 1 : 10
%     total_snr15{i} = [];
% end
% 
% % Initialize a counter
% counter = 0;
% 
% % Outter for loop
% for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
%         'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
%         'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
%         'SX397.WAV'}
%     % Read the PHN TIMIT file
%     [s,fs,wrd,phn] = readsph(filename2{1},'wt');
% 
%     % Use counter 
%     counter = counter + 1;
%     
%     % Inner for loop
%     for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
%             'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
%             'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
%             'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
%         % Read the noise file
%         [s2, fs2, bits] = wavread(filename{1});
% 
%         % Define the SNR
%         snr = 15;
% 
%         % Add noise from column vector n() at sample frequency fn with random start
%         % Sample and wrapping around as required with a vorbis cross-fade
%         z = v_addnoise(s,fs,snr,'',s2,fs2); 
% 
%         % store the result
%         total_snr15{counter} = [total_snr15{counter} z];
%     end
% end
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Define the matrix for SNR 10
% for i = 1 : 10
%     total_snr10{i} = [];
% end
% 
% % Initialize a counter
% counter = 0;
% 
% % Outter for loop
% for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
%         'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
%         'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
%         'SX397.WAV'}
%     % Read the PHN TIMIT file
%     [s,fs,wrd,phn] = readsph(filename2{1},'wt');
% 
%     % Use counter 
%     counter = counter + 1;
%     
%     % Inner for loop
%     for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
%             'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
%             'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
%             'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
%         % Read the noise file
%         [s2, fs2, bits] = wavread(filename{1});
% 
%         % Define the SNR
%         snr = 10;
% 
%         % Add noise from column vector n() at sample frequency fn with random start
%         % Sample and wrapping around as required with a vorbis cross-fade
%         z = v_addnoise(s,fs,snr,'',s2,fs2); 
% 
%         % store the result
%         total_snr10{counter} = [total_snr10{counter} z];
%     end
% end
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Define the matrix for SNR 5
% for i = 1 : 10
%     total_snr5{i} = [];
% end
% 
% % Initialize a counter
% counter = 0;
% 
% % Outter for loop
% for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
%         'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
%         'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
%         'SX397.WAV'}
%     % Read the PHN TIMIT file
%     [s,fs,wrd,phn] = readsph(filename2{1},'wt');
% 
%     % Use counter 
%     counter = counter + 1;
%     
%     % Inner for loop
%     for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
%             'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
%             'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
%             'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
%         % Read the noise file
%         [s2, fs2, bits] = wavread(filename{1});
% 
%         % Define the SNR
%         snr = 10;
% 
%         % Add noise from column vector n() at sample frequency fn with random start
%         % Sample and wrapping around as required with a vorbis cross-fade
%         z = v_addnoise(s,fs,snr,'',s2,fs2); 
% 
%         % store the result
%         total_snr5{counter} = [total_snr5{counter} z];
%     end
% end
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Define the matrix for SNR 0
% for i = 1 : 10
%     total_snr0{i} = [];
% end
% 
% % Initialize a counter
% counter = 0;
% 
% % Outter for loop
% for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
%         'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
%         'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
%         'SX397.WAV'}
%     % Read the PHN TIMIT file
%     [s,fs,wrd,phn] = readsph(filename2{1},'wt');
% 
%     % Use counter 
%     counter = counter + 1;
%     
%     % Inner for loop
%     for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
%             'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
%             'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
%             'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
%         % Read the noise file
%         [s2, fs2, bits] = wavread(filename{1});
% 
%         % Define the SNR
%         snr = 10;
% 
%         % Add noise from column vector n() at sample frequency fn with random start
%         % Sample and wrapping around as required with a vorbis cross-fade
%         z = v_addnoise(s,fs,snr,'',s2,fs2); 
% 
%         % store the result
%         total_snr0{counter} = [total_snr0{counter} z];
%     end
% end
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Store the variables
% save('noisy_speech_TIMIT_TRAIN_FCJF0_snr15.mat', 'total_snr15');
% save('noisy_speech_TIMIT_TRAIN_FCJF0_snr10.mat', 'total_snr10');
% save('noisy_speech_TIMIT_TRAIN_FCJF0_snr5.mat', 'total_snr5');
% save('noisy_speech_TIMIT_TRAIN_FCJF0_snr0.mat', 'total_snr0');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% The noisy signals have been created

% We load the signals
load('noisy_speech_TIMIT_TRAIN_FCJF0_snr15.mat');
load('noisy_speech_TIMIT_TRAIN_FCJF0_snr10.mat');
load('noisy_speech_TIMIT_TRAIN_FCJF0_snr5.mat');
load('noisy_speech_TIMIT_TRAIN_FCJF0_snr0.mat');

% We can use: "soundsc(total_snr15{1}(:,1), 16000)"
% We can use: "soundsc(total_snr10{1}(:,1), 16000)"
% We can use: "soundsc(total_snr5{1}(:,1), 16000)"
% We can use: "soundsc(total_snr0{1}(:,1), 16000)"
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % Define the matrix for SNR -5
% for i = 1 : 10
%     total_snr_minus5{i} = [];
% end
% 
% % Initialize a counter
% counter = 0;
% 
% % Outter for loop
% for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
%         'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
%         'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
%         'SX397.WAV'}
%     % Read the PHN TIMIT file
%     [s,fs,wrd,phn] = readsph(filename2{1},'wt');
% 
%     % Use counter 
%     counter = counter + 1;
%     
%     % Inner for loop
%     for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
%             'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
%             'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
%             'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
%         % Read the noise file
%         [s2, fs2, bits] = wavread(filename{1});
% 
%         % Define the SNR
%         snr = -5;
% 
%         % Add noise from column vector n() at sample frequency fn with random start
%         % Sample and wrapping around as required with a vorbis cross-fade
%         z = v_addnoise(s,fs,snr,'',s2,fs2); 
% 
%         % store the result
%         total_snr_minus5{counter} = [total_snr_minus5{counter} z];
%     end
% end
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Define the matrix for SNR -10
% for i = 1 : 10
%     total_snr_minus10{i} = [];
% end
% 
% % Initialize a counter
% counter = 0;
% 
% % Outter for loop
% for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
%         'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
%         'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
%         'SX397.WAV'}
%     % Read the PHN TIMIT file
%     [s,fs,wrd,phn] = readsph(filename2{1},'wt');
% 
%     % Use counter 
%     counter = counter + 1;
%     
%     % Inner for loop
%     for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
%             'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
%             'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
%             'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
%         % Read the noise file
%         [s2, fs2, bits] = wavread(filename{1});
% 
%         % Define the SNR
%         snr = -10;
% 
%         % Add noise from column vector n() at sample frequency fn with random start
%         % Sample and wrapping around as required with a vorbis cross-fade
%         z = v_addnoise(s,fs,snr,'',s2,fs2); 
% 
%         % store the result
%         total_snr_minus10{counter} = [total_snr_minus10{counter} z];
%     end
% end
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Store the variables
% save('noisy_speech_TIMIT_TRAIN_FCJF0_snr_minus5.mat', 'total_snr_minus5');
% save('noisy_speech_TIMIT_TRAIN_FCJF0_snr_minus10.mat', 'total_snr_minus10');
% % We use: "soundsc(total_snr_minus10{1}(:,1),16000)"
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% The noisy signals have been created

% We load the signals
load('noisy_speech_TIMIT_TRAIN_FCJF0_snr_minus5.mat');
load('noisy_speech_TIMIT_TRAIN_FCJF0_snr_minus10.mat');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We now find a model for high noise power speech

% % We use: "total_snr_minus10", "total_snr_minus5"
% % We use: "total_snr5", "total_snr0"
% 
% % We use: "total_snr_minus10{1}(:,1)"
% % We use: "total_snr_minus5{1}(:,1)"
% % We can use: "total_snr5{1}(:,1)"
% % We can use: "total_snr0{1}(:,1)"
% 
% % We initialize the matrix.
% total_speech_active_frames_noisy = [];
% 
% % Outter for loop
% for outter_loop = 1 : 10
%     for inner_loop = 1 : 15
%        
%         %[x, Srate, bits] = wavread( filename);	
%         x = total_snr_minus10{outter_loop}(:, inner_loop);
% 
%         % Define Fs
%         Srate = 16000;
%         
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%             % Speech and noise frame found
%             else
%                 total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
%         %end
%         end
%     end
% end
% 
% % Store the variable
% save('frames_total_snr_minus10.mat', 'total_speech_active_frames_noisy');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We load the voiced frames
frames1 = load('frames_total_snr_minus10.mat');
frames2 = load('frames_total_snr_minus5.mat');
frames3 = load('frames_total_snr5.mat');
frames4 = load('frames_total_snr0.mat');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use: "frames1.total_speech_active_frames_noisy"
% We use: "frames2.total_speech_active_frames_noisy"
% We use: "frames3.total_speech_active_frames_noisy"
% We use: "frames4.total_speech_active_frames_noisy"

% Initialize the matrix with the noisy speech frames
final_noisy_frames = [];

% Define the matrix with the noisy speech frames
final_noisy_frames = [frames1.total_speech_active_frames_noisy; ...
    frames2.total_speech_active_frames_noisy; ...
    frames3.total_speech_active_frames_noisy; ...
    frames4.total_speech_active_frames_noisy];

% Define the frame size in samples
len = size(final_noisy_frames,2); 

% Define the window
win = hanning(len)';  

% normalize window for equal level output 
win = win*len/sum(win);  

% Define the matrix for the power at 800 Hz
power_at_800Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define N
N = size(final_noisy_frames,2);

% for loop for all frames
for i = 1 : size(final_noisy_frames,1)
    % We use the first frame: "voiced_clean_speech_frames.total_speech_active_frames(i,:)"
    insign = win .* final_noisy_frames(i,:);
    % We have used windowing.

    % Perform the FFT
    xdft = fft(insign);
    xdft = xdft(1:N/2+1);

    % Find the power spectral density
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    % Define the frequency axis
    freq = 0:Fs/length(insign):Fs/2;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the speech level in units of power
    lev = activlev(final_noisy_frames(i,:), Fs);     
    
    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)/lev];
end

% plot the histogram
nbins = 100;
[y,x] = hist(log10(power_at_800Hz_noisy), nbins);
%[y,x] = hist(power_at_800Hz_noisy, nbins);

% Plot the histogram
figure; set(gcf,'Color','w');
bar(x, y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the voiced (very) noisy speech power at 800 Hz. Log power is used.');
ylabel('Counts'); xlabel('Log power: 10 * log10(power)');
%ylabel('Counts'); xlabel('Power');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use GMM
data = 10*log10(power_at_800Hz_noisy);

% create GMM with k mixtures
x = data'; k = 2;
[m,v,w] = gaussmix2(x,[],[],k,'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 3;
[m2,v2,w2] = gaussmix2(x,[],[],k, 'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 4;
[m3,v3,w3] = gaussmix2(x,[],[],k, 'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% We model high power noise 
% We use noise only frames

% Load files
load('noisy_speech_TIMIT_TRAIN_FCJF0_snr_minus5.mat');
load('noisy_speech_TIMIT_TRAIN_FCJF0_snr_minus10.mat');

% % We use: "total_snr_minus10", "total_snr_minus5"
% % We use: "total_snr5", "total_snr0"
% 
% % We use: "total_snr_minus10{1}(:,1)"
% % We use: "total_snr_minus5{1}(:,1)"
% % We can use: "total_snr5{1}(:,1)"
% % We can use: "total_snr0{1}(:,1)"
% 
% % We initialize the matrix.
% total_speech_active_frames_noisy = [];
% 
% % Outter for loop
% for outter_loop = 1 : 10
%     for inner_loop = 1 : 15
%        
%         %[x, Srate, bits] = wavread( filename);	
%         x = total_snr5{outter_loop}(:, inner_loop);
% 
%         % Define Fs
%         Srate = 16000;
%         
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%                 total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
%         %end
%         end
%     end
% end
% 
% % Store the variable
% save('frames_noiseOnly_total_snr5.mat', 'total_speech_active_frames_noisy');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Load data
noiseOnly1 = load('frames_noiseOnly_total_snr_minus10.mat');
noiseOnly2 = load('frames_noiseOnly_total_snr_minus5.mat');
noiseOnly3 = load('frames_noiseOnly_total_snr5.mat');
noiseOnly4 = load('frames_noiseOnly_total_snr0.mat');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Initialize the matrix with the noisy speech frames
final_noisy_frames = [];

% Define the matrix with the noisy speech frames
 final_noisy_frames = [noiseOnly1.total_speech_active_frames_noisy; ...
     noiseOnly2.total_speech_active_frames_noisy; ...
     noiseOnly3.total_speech_active_frames_noisy; ...
     noiseOnly4.total_speech_active_frames_noisy];
%final_noisy_frames = noiseOnly4.total_speech_active_frames_noisy;

% Define the frame size in samples
len = size(final_noisy_frames,2); 

% Define the window
win = hanning(len)';  

% normalize window for equal level output 
win = win*len/sum(win);  

% Define the matrix for the power at 800 Hz
power_at_800Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define N
N = size(final_noisy_frames,2);

% for loop for all frames
for i = 1 : size(final_noisy_frames,1)
    % We use the first frame: "voiced_clean_speech_frames.total_speech_active_frames(i,:)"
    insign = win .* final_noisy_frames(i,:);
    % We have used windowing.

    % Perform the FFT
    xdft = fft(insign);
    xdft = xdft(1:N/2+1);

    % Find the power spectral density
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    % Define the frequency axis
    freq = 0:Fs/length(insign):Fs/2;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the speech level in units of power
    lev = activlev(final_noisy_frames(i,:), Fs);     
    
    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)/lev];
end

% plot the histogram
nbins = 100;
%[y,x] = hist(log10(power_at_800Hz_noisy), nbins);
[y,x] = hist(power_at_800Hz_noisy, nbins);

% Plot the histogram
figure; set(gcf,'Color','w');
bar(x, y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the high noise power at 800 Hz. Log power is not used.');
%ylabel('Counts'); xlabel('Log power: 10 * log10(power)');
ylabel('Counts'); xlabel('Power');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use GMM
data = power_at_800Hz_noisy;

% create GMM with k mixtures
x = data'; k = 2;
[m,v,w] = gaussmix2(x,[],[],k,'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 3;
[m2,v2,w2] = gaussmix2(x,[],[],k, 'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 4;
[m3,v3,w3] = gaussmix2(x,[],[],k, 'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% iteration, loop for the 10 different model orders
data = power_at_800Hz_noisy;
loop_AIC = 0; 
N = length(data);

%we find the AIC and MDL for each case using the following code
for q = 1 : 30 
    % We use GMM
    % Create an object with this GMM
    GMFIT = gmdistribution.fit(data', q);

    % We use the AIC
    loop_AIC(q,1) = GMFIT.AIC;
end

% plot the AIC with q, q = [1:10]
figure; set(gcf,'Color','w');
plot([1:30], loop_AIC);
grid on; grid minor;

% name the axes
set(gca,'FontSize',15);
title('AIC for GMM'); xlabel('Model Order p '); ylabel('AIC');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use k = 6
k = 6;
GMFIT = gmdistribution.fit(data', k);
%[m, v, w] = gaussmix(data', [], [], k);    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% We model low power noise
 
% % We use: "total_snr15", "total_snr10"
% % We use: "total_snr5"
% 
% % We use: "total_snr15{1}(:,1)"
% % We can use: "total_snr10{1}(:,1)"
% % We can use: "total_snr5{1}(:,1)"
% 
% % We initialize the matrix.
% total_speech_active_frames_noisy = [];
% 
% % Outter for loop
% for outter_loop = 1 : 10
%     for inner_loop = 1 : 15
%        
%         %[x, Srate, bits] = wavread( filename);	
%         x = total_snr5{outter_loop}(:, inner_loop);
% 
%         % Define Fs
%         Srate = 16000;
%         
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%                 total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
%         %end
%         end
%     end
% end
% 
% % Store the variable
% save('frames2_noiseOnly_total_snr5.mat', 'total_speech_active_frames_noisy');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Load data
new_noiseOnly1 = load('frames2_noiseOnly_total_snr5.mat');
new_noiseOnly2 = load('frames_noiseOnly_total_snr10.mat');
new_noiseOnly3 = load('frames_noiseOnly_total_snr15.mat');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Initialize the matrix with the noisy speech frames
final_noisy_frames = [];

% Define the matrix with the noisy speech frames
final_noisy_frames = [new_noiseOnly1.total_speech_active_frames_noisy; ...
     new_noiseOnly2.total_speech_active_frames_noisy; ...
     new_noiseOnly3.total_speech_active_frames_noisy];
%final_noisy_frames = new_noiseOnly3.total_speech_active_frames_noisy;

% Define the frame size in samples
len = size(final_noisy_frames,2); 

% Define the window
win = hanning(len)';  

% normalize window for equal level output 
win = win*len/sum(win);  

% Define the matrix for the power at 800 Hz
power_at_800Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define N
N = size(final_noisy_frames,2);

% for loop for all frames
for i = 1 : size(final_noisy_frames,1)
    % We use the first frame: "voiced_clean_speech_frames.total_speech_active_frames(i,:)"
    insign = win .* final_noisy_frames(i,:);
    % We have used windowing.

    % Perform the FFT
    xdft = fft(insign);
    xdft = xdft(1:N/2+1);

    % Find the power spectral density
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    % Define the frequency axis
    freq = 0:Fs/length(insign):Fs/2;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the speech level in units of power
    lev = activlev(final_noisy_frames(i,:), Fs);     
    
    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)/lev];
end

% plot the histogram
nbins = 100;
%[y,x] = hist(log10(power_at_800Hz_noisy), nbins);
[y,x] = hist(power_at_800Hz_noisy, nbins);

% Plot the histogram
figure; set(gcf,'Color','w');
bar(x, y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the low noise power at 800 Hz. Log power is not used.');
ylabel('Counts'); xlabel('Power');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use GMM
%data = 10*log10(power_at_800Hz_noisy);
data = power_at_800Hz_noisy;

% create GMM with k mixtures
x = data'; k = 2;
[m,v,w] = gaussmix2(x,[],[],k,'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 3;
[m2,v2,w2] = gaussmix2(x,[],[],k, 'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 4;
[m3,v3,w3] = gaussmix2(x,[],[],k, 'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% iteration, loop for the 10 different model orders
data = power_at_800Hz_noisy;
loop_AIC = 0; 
N = length(data);

%we find the AIC and MDL for each case using the following code
for q = 1 : 30 
    % We use GMM
    % Create an object with this GMM
    GMFIT = gmdistribution.fit(data', q);

    % We use the AIC
    loop_AIC(q,1) = GMFIT.AIC;
end

% plot the AIC with q, q = [1:10]
figure; set(gcf,'Color','w');
plot([1:30], loop_AIC);
grid on; grid minor;

% name the axes
set(gca,'FontSize',15);
title('AIC for GMM'); xlabel('Model Order p '); ylabel('AIC');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use k = 6
k = 6;
GMFIT = gmdistribution.fit(data', k);
%[m, v, w] = gaussmix(data', [], [], k);    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------













%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% We use the PEFAC algorithm

% Read the PHN TIMIT file
filename2 = 'SA1.WAV';
% [s,fs,wrd,phn] = readsph(filename2,'wt');

[s,fs] = readsph((filename2));

% Use "wavwrite"
%wavwrite(y, fs, bits, 'asdfasdfa');
%writewav(s,fs, '1asfdasfdasa');

% Read the noise file
filename = 'babble'; 
[s2,fs2] = readwav((filename));
%f = v_resample(s2,4000,4995);

% Resample
%[y,~] = v_resample(s2, 999,800);

% Use "wavwrite"
%wavwrite(y, fs, bits, 'asdfasdfa');
%writewav(y,fs, 'asdfasdfa');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'AD',s2,fs2); 

% store z
% Store the result to a file 
filename_noisy = 'noisy_speech_signal_file2.wav';
wavwrite(z, fs, filename_noisy);

% or:  z=v_addnoise(s,fs,snr,'AD',s2,fs2);

%fxpefac(z,16000)
figure; set(gcf,'Color','w');
fxpefac(z, fs);
set(gca,'FontSize',15);
title('Without pre-processing. Results of the PEFAC algorithm. SNR = 5 dB.');

% new
% we use "z"

% we subtract transient noise
%xfinal = SS_transient_noise( z );
% test the new MS algorithm

% Read the noisy speech file
filename = filename_noisy;
% [input_main_noisy_speech_signal_z, ~] = readwav((filename));

% (5) Minimum tracking algorithm, Martin (2001), (use: 'martin')
% Define the method to be used 
method = 'martin';

% Call the function so as to find the noise power spectrum for each frame
[martin_noise_ps, minact_new] = new_noise_parameters(filename, method);

% implement the new algorithm

% (5) Minimum tracking algorithm, Martin (2001), (use: 'martin')
% Define the method to be used 
method = 'martin2';

% Call the function so as to find the noise power spectrum for each frame
[final_clean_speech, final_noise_only] = new_SS_specsub_ns(filename, method, minact_new);
xfinal = final_clean_speech;

%fxpefac(z,16000)
figure; set(gcf,'Color','w');
fxpefac(xfinal, fs);
set(gca,'FontSize',15);
title('With pre-processing. Results of the PEFAC algorithm. SNR = 5 dB.');

%fxpefac(z,16000)
figure; set(gcf,'Color','w');
fxpefac(s, fs);
set(gca,'FontSize',15);
title('True PEFAC. No noise. Results of the PEFAC algorithm. SNR = Infinity.');

% Define the SNR
snr = 0;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'AD',s2,fs2); 

% or:  z=v_addnoise(s,fs,snr,'AD',s2,fs2);

%fxpefac(z,16000)
figure; set(gcf,'Color','w');
fxpefac(z, fs);
set(gca,'FontSize',15);
title('Without pre-processing. Results of the PEFAC algorithm. SNR = 0 dB.');

% find the active level of speech
activlevg(z, fs)
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % read the file
% filename = 'sp01_16kHz';
% %[s, Fs, ~] = wavread(filename);
% [s, Fs, ~] = readwav(filename);
% 
% % Plot the PEFAC result
% figure; set(gcf,'Color','w');
% fxpefac(s, Fs);
% 
% % Set the title of the image
% set(gca,'FontSize',15);
% title('Results of the PEFC algorithm. No noise is used.');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % Read the noise file
% filename = 'sp01_babble_sn0'; 
% %[s2, fs2, bits] = wavread(filename);
% [s2, fs2, bits] = readwav(filename);
% 
% % Define the SNR
% snr = 10;
% 
% % Add noise from column vector n() at sample frequency fn with random start
% % Sample and wrapping around as required with a vorbis cross-fade
% z = v_addnoise(s,Fs,snr,'',s2,fs2); 
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
 
% % read the file
% filename = 'sp01_airport_sn15';
% [s, Fs, ~] = wavread(filename);
% 
% % Plot the PEFAC result
% figure; set(gcf,'Color','w');
% fxpefac(z, Fs);
% 
% % Set the title of the image
% set(gca,'FontSize',15);
% title('Results of the PEFC algorithm. SNR 0 is used.');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
 
% % Plot the PEFAC result
% figure; set(gcf,'Color','w');
% fxpefac(total_snr0{1}(:,1), 16000);
% 
% % Set the title of the image
% set(gca,'FontSize',15);
% title('Results of the PEFC algorithm. SNR 0 is used.');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
 
% % We use the PEFAC algortihm
% s = total_snr15{1}(:,1);
% fs = 16000;
% 
% figure; set(gcf,'Color','w');
% fxpefac(s,fs);
% 
% figure; set(gcf,'Color','w');
% [fx,tx,pv,fv] = fxpefac(s,fs, 0.01, 'g');
% 
% 
% [fx,tx,pv,fv] = fxpefac(s,fs);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
 
% filename = 'sp01_16kHz';
% [s, sdasd, sa] = wavread(filename);
% 
% sp = s;
% fs= 16000;
% lev = activlevg(sp,fs, 'n');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------


















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%

% this is method 1
% we use the Normal Distribution

% Read the PHN TIMIT file
filename2 = 'SA2.WAV';
% [s,fs,wrd,phn] = readsph(filename2,'wt');

[s,fs] = readsph((filename2));

% Use "wavwrite"
%wavwrite(y, fs, bits, 'asdfasdfa');
%writewav(s,fs, '1asfdasfdasa');

% Read the noise file
filename = 'machinegun'; 
[s2,fs2] = readwav((filename));
%f = v_resample(s2,4000,4995);

% Resample
%[y,~] = v_resample(s2, 999,800);

% Use "wavwrite"
%wavwrite(y, fs, bits, 'asdfasdfa');
%writewav(y,fs, 'asdfasdfa');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

% We use "z", we use "fs"
% Define Srate
Srate = fs;
% Define signal
x = z;

% Define the frame size in samples
len=floor(20*Srate/1000); 
if rem(len,2)==1, len=len+1; end;

% window overlap in percent of frame size
PERC=50; 
len1=floor(len*PERC/100);
len2=len-len1;

% define window
win=hanning(len);  

% normalize window for equal level output 
win = win*len2/sum(win);  

% Noise magnitude calculations - assuming that the first 6 frames is noise/silence
nFFT=2*len;
j=1;
noise_mean=zeros(nFFT,1);
for k=1:6
    noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
    j=j+len;
end
noise_mu=noise_mean/6;
noise_mu2=noise_mu.^2;

%--- allocate memory and initialize various variables
k=1;
img=sqrt(-1);
x_old=zeros(len1,1);
Nframes=floor(length(x)/len2)-1;
xfinal=zeros(Nframes*len2,1);

% --------------- Initialize parameters ------------
k=1;
aa=0.98;
eta= 0.15;
mu=0.98;
c=sqrt(pi)/2;
qk=0.3;
qkr=(1-qk)/qk;
ksi_min=10^(-25/10); 

power_noise = [];
index_power = [];

% Main for loop
for n = 1 : Nframes

    insign=win.*x(k:k+len-1);

    %--- Take fourier transform of  frame
    spec=fft(insign,nFFT);
    sig=abs(spec); % compute the magnitude
    sig2=sig.^2;

    gammak=min(sig2./noise_mu2,40);  % posteriori SNR
    if n==1
        ksi=aa+(1-aa)*max(gammak-1,0);
    else
        ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);   
        
        % decision-direct estimate of a priori SNR
        ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
    end

    log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
    vad_decision = sum(log_sigma_k)/nFFT;    

    % Noise only frame found
    if (vad_decision < eta) 
        % store power 
        power_noise = [power_noise; spec'];
        
        % store the index
        index_power = [index_power n];

        % noise_mu2 = mu * noise_mu2 + (1 - mu) * sig2;
                
%     % Speech and noise frame found
%     else
%         total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
    end

    vk=ksi.*gammak./(1+ksi);
    j0 = besseli(0,vk/2);
    j1 = besseli(1,vk/2);

    C=exp(-0.5*vk);
    A=((c*(vk.^0.5)).*C)./gammak;
    B=(1+vk).*j0+vk.*j1;
    hw=A.*B;

    sig=sig.*hw;

    % save for estimation of a priori SNR in next frame
    Xk_prev=sig.^2;  

    xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);

    xi_w= real( xi_w);

    xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
    x_old= xi_w(len1+ 1: len);

    k=k+len2;
%end
end

[a, b] = size(power_noise);
variance_hat = [];

power_noise_real = real(power_noise);

for i = 1 : a
    % the variance is the same in the real and imaginary parts    
    data = power_noise_real(i,:);
    [~, sigma_hat] = normfit(data);
    
    % store the variance
    variance_hat = [variance_hat sigma_hat^2];
end

% the variance is the same in the real and imaginary parts    

% we find the mean power parameter
% the mean power parameter is the only parameter of the Negative Exponential Distribution
variance_hat = 2 * variance_hat;

% we find the mean power parameter
mean_power_parameter = variance_hat;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% this is method 2
% we use 3 parameters, we use a GMM with 2 Gaussians

% Read the PHN TIMIT file
filename2 = 'SA2.WAV';
% [s,fs,wrd,phn] = readsph(filename2,'wt');

[s,fs] = readsph((filename2));

% Use "wavwrite"
%wavwrite(y, fs, bits, 'asdfasdfa');
%writewav(s,fs, '1asfdasfdasa');

% Read the noise file
filename = 'machinegun'; 
[s2,fs2] = readwav((filename));
%f = v_resample(s2,4000,4995);

% Resample
%[y,~] = v_resample(s2, 999,800);

% Use "wavwrite"
%wavwrite(y, fs, bits, 'asdfasdfa');
%writewav(y,fs, 'asdfasdfa');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

% We use "z", we use "fs"
% Define Srate
Srate = fs;
% Define signal
x = z;

% Define the frame size in samples
len=floor(20*Srate/1000); 
if rem(len,2)==1, len=len+1; end;

% window overlap in percent of frame size
PERC=50; 
len1=floor(len*PERC/100);
len2=len-len1;

% define window
win=hanning(len);  

% normalize window for equal level output 
win = win*len2/sum(win);  

% Noise magnitude calculations - assuming that the first 6 frames is noise/silence
nFFT=2*len;
j=1;
noise_mean=zeros(nFFT,1);
for k=1:6
    noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
    j=j+len;
end
noise_mu=noise_mean/6;
noise_mu2=noise_mu.^2;

%--- allocate memory and initialize various variables
k=1;
img=sqrt(-1);
x_old=zeros(len1,1);
Nframes=floor(length(x)/len2)-1;
xfinal=zeros(Nframes*len2,1);

% --------------- Initialize parameters ------------
k=1;
aa=0.98;
eta= 0.15;
mu=0.98;
c=sqrt(pi)/2;
qk=0.3;
qkr=(1-qk)/qk;
ksi_min=10^(-25/10); 

power_noise = [];
index_power2 = [];

% Main for loop
for n = 1 : Nframes

    insign=win.*x(k:k+len-1);

    %--- Take fourier transform of  frame
    spec=fft(insign,nFFT);
    sig=abs(spec); % compute the magnitude
    sig2=sig.^2;

    gammak=min(sig2./noise_mu2,40);  % posteriori SNR
    if n==1
        ksi=aa+(1-aa)*max(gammak-1,0);
    else
        ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);   
        
        % decision-direct estimate of a priori SNR
        ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
    end

    log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
    vad_decision = sum(log_sigma_k)/nFFT;    

    % Noise only frame found
    if (vad_decision < eta) 
        % noise_mu2 = mu * noise_mu2 + (1 - mu) * sig2;
        
        % store power 
        power_noise = [power_noise; spec'];
        
        % store the index
        index_power2 = [index_power2 n];
        
%     % Speech and noise frame found
%     else
%         total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
    end

    vk=ksi.*gammak./(1+ksi);
    j0 = besseli(0,vk/2);
    j1 = besseli(1,vk/2);

    C=exp(-0.5*vk);
    A=((c*(vk.^0.5)).*C)./gammak;
    B=(1+vk).*j0+vk.*j1;
    hw=A.*B;

    sig=sig.*hw;

    % save for estimation of a priori SNR in next frame
    Xk_prev=sig.^2;  

    xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);

    xi_w= real( xi_w);

    xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
    x_old= xi_w(len1+ 1: len);

    k=k+len2;
%end
end

[a, b] = size(power_noise);
variance_hat = [];
weight_hat = [];

power_noise_real = real(power_noise);

for i = 1 : a
    % the variance is the same in the real and imaginary parts    
    data = power_noise_real(i,:);
    
    % create GMM with k = 2 mixtures
    x = data'; 
    k = 2;
    [~, v, w] = gaussmix(x,[],[],k, 'vkf');
    % [~, v, w] = gaussmix(x,[],[],k,'vf');   
    
    % store the variance
    variance_hat = [variance_hat; [v(:,:,1) v(:,:,2)]];

    % store the weights
    weight_hat = [weight_hat; w'];
end

% we delete the second row of the weights
% the weights sum to 1
weight_hat = weight_hat(:,1);

% the variance is the same in the real and imaginary parts    
% we find the mean power parameter
variance_hat = 2 * variance_hat;

% the 3 parameters are: variance_hat, weight_hat

% we find the mean power parameter
mean_power_parameter2 = [];
for i = 1 : length(weight_hat)
    mean_power_parameter2(i,1) = (weight_hat(i) * variance_hat(i,1)) + ((1-weight_hat(i)) * variance_hat(i,2));
end
mean_power_parameter2 = mean_power_parameter2';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% VAD: we update the noise power only when no speech 

% define 2 matrices
final_noise_power = zeros(1,Nframes);
final_noise_power2 = zeros(1,Nframes);

for i = 1 : Nframes
    for k = 1 : length(index_power)
        if (i == index_power(k))
            final_noise_power(i) = mean_power_parameter(k);
        end
    end
end

for i = 1 : Nframes
    if (final_noise_power(i) == 0)
        final_noise_power(i) = final_noise_power(i-1);
    end
end


for i = 1 : Nframes
    for k = 1 : length(index_power2)
        if (i == index_power2(k))
            final_noise_power2(i) = mean_power_parameter2(k);
        end
    end
end

for i = 1 : Nframes
    if (final_noise_power2(i) == 0)
        final_noise_power2(i) = final_noise_power2(i-1);
    end
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------





















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% We now use the "NOIZEUS_Noisy_Speech_Signal_Database_16kHz" database.
% We will use the noisy speech signals

% % Define the file name
% filename = 'sp01_train_sn0';
% 
% % We initialize the matrix.
% total_speech_frames_noisy = [];
% 
% % Outter for loop
% for outter_loop = 1 : 4
%     for inner_loop = 1 : 30
%         % We change the filename based on inner_loop
%         % We use a temporary variable
%         temp_variable = num2str(inner_loop);
% 
%         % We now change the filename based on inner_loop
%         if (inner_loop < 10) 
%             filename(3) = num2str(0);
%             filename(4) = num2str(temp_variable(1));
%         else
%             filename(3) = num2str(temp_variable(1));
%             filename(4) = num2str(temp_variable(2));
%         end
% 
%         % If statement to change the filename
%         if (outter_loop == 1) 
%             filename(14) = num2str(0);
%         elseif (outter_loop == 2) 
%             filename(14) = num2str(5);
%         elseif (outter_loop == 3) 
%             filename(14) = num2str(1);
%             filename(15) = num2str(0);
%         else 
%             filename(14) = num2str(1);
%             filename(15) = num2str(5);
%         end
% 
%         %if (strcmp(filename, 'sp04_babble_sn10') ~= 1)
%         % Read the .wav file
%         [x, Srate, bits] = wavread( filename);	
% 
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             % Store the power spectrum
%             % The variable "sig2" is of size 640 x 1
%             total_speech_frames_noisy = [total_speech_frames_noisy sig2];
%         %end
%         end
%     end
% end
% 
% % Store the variable
% save('noise_speech_8.mat', 'total_speech_frames_noisy');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Load the data
noise1 = load('noise_speech_1.mat');
noise2 = load('noise_speech_2.mat');
noise3 = load('noise_speech_3.mat');
noise4 = load('noise_speech_4.mat');
noise5 = load('noise_speech_5.mat');
noise6 = load('noise_speech_6.mat');
noise7 = load('noise_speech_7.mat');
noise8 = load('noise_speech_8.mat');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

power_at_800Hz_noisy = [];
Fs = 16000;

total_speech_frames_noisy = [noise1.total_speech_frames_noisy, noise2.total_speech_frames_noisy, ...
    noise3.total_speech_frames_noisy, noise4.total_speech_frames_noisy, noise5.total_speech_frames_noisy, ...
    noise6.total_speech_frames_noisy, noise7.total_speech_frames_noisy, noise8.total_speech_frames_noisy];

for i = 1 : size(total_speech_frames_noisy ,2)
    psdx = total_speech_frames_noisy(:,i);

    freq = 0:Fs/(640*2):Fs/2-1;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)];
end

% plot the histogram
nbins = 100;
[y,x] = hist(10*log10(power_at_800Hz_noisy), nbins);

% Plot the histogram
figure; set(gcf,'Color','w');
bar(x, y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the noisy speech power at 800 Hz. Log power is used.');
ylabel('Counts'); xlabel('Log power: 10 * log10(power)');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------














%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% We use a VAD to find the clean speech voiced frames.

% % Define the filename
% filename = 'zpg_in_1_1';
% 
% % Initialize the variables for the for loop
% total_speech_active_frames = [];
% 
% % Outter for loop
% for outter_loop = 1 : 9
%     % Inner for loop
%     for inner_loop = 1 : 3
%         [x, Fs, nbits] = wavread(filename);
% 
%         % We change the filename based on inner_loop
%         filename(4+6) = num2str(inner_loop);
% 
%         % We change the filename based on outter_loop
%         filename(4+4) = num2str(outter_loop);
% 
%         % Read the .wav file
%         [x, Srate, bits] = wavread( filename);	
% 
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%             % Speech and noise frame found
%             else
%                 total_speech_active_frames = [total_speech_active_frames; (x(k:k+len-1))'];
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
% 
%         end
%     end
% end
% 
% % Define the filename
% filename = 'zpg_in_10_1';
% 
% % Outter for loop
% for outter_loop = 10 : 50
%     % Inner for loop
%     for inner_loop = 1 : 3
%         [x, Fs, nbits] = wavread(filename);
% 
%         % We change the filename based on inner_loop
%         filename(4+7) = num2str(inner_loop);
% 
%         % We use a temporary variable
%         temp_variable = num2str(outter_loop);
%     
%         % We now change the filename based on k
%         filename(4+4) = num2str(temp_variable(1));
%         filename(4+5) = num2str(temp_variable(2));
% 
%         % Read the .wav file
%         [x, Srate, bits] = wavread(filename);	
% 
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%             % Speech and noise frame found
%             else
%                 total_speech_active_frames = [total_speech_active_frames; (x(k:k+len-1))'];
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
% 
%         end
%     end
% end
% 
% % Store the variable
% save('clean_speech_frames_zpg.mat', 'total_speech_active_frames');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % Read the .wav file
% filename = 'jal_in_1_1';
% [x, Srate, bits] = wavread(filename);	
% 
% % Play the audio file
% %wavplay(x,Srate);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Load all the data

% Initialize the variable
total_speech_frames_noisy = [];

%clean1 = load('clean_speech_frames.mat');
clean2 = load('clean_speech_frames_jal.mat');
total_speech_frames_noisy = clean2.total_speech_active_frames';

clean2 = load('clean_speech_frames_jcs.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];

clean2 = load('clean_speech_frames_jgl.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_leb.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_mam.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_mwm.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_nad.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_sas.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_scs.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_sll.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_zng.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
clean2 = load('clean_speech_frames_zpg.mat');
total_speech_frames_noisy = [total_speech_frames_noisy clean2.total_speech_active_frames'];
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

power_at_800Hz_noisy = [];
Fs = 16000;

%total_speech_frames_noisy = [clean2.total_speech_active_frames', ...
%    clean3.total_speech_active_frames', clean4.total_speech_active_frames', clean5.total_speech_active_frames', ...
%    clean6.total_speech_active_frames', clean7.total_speech_active_frames', clean8.total_speech_active_frames'];, ...
%    clean9.total_speech_active_frames', clean10.total_speech_active_frames', clean11.total_speech_active_frames', ...
%    clean12.total_speech_active_frames', clean13.total_speech_active_frames'];

for i = 1 : size(total_speech_frames_noisy ,2)
    psdx = total_speech_frames_noisy(:,i);

    freq = 0:Fs/(640*2):Fs/2-1;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)];
end

% plot the histogram
figure; set(gcf,'Color','w');
%nhist(10*log10(power_at_800Hz_noisy));
nhist(power_at_800Hz_noisy);

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the noisy speech power at 800 Hz. Log power is not used.');
ylabel('Counts'); xlabel('Power');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

nbins = 100;
[y,x] = hist(power_at_800Hz_noisy, nbins);
% Plot the histogram
bar(x, y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the noisy speech power at 800 Hz. Log power is not used.');
ylabel('Counts'); xlabel('Power');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

data = 10*log10(power_at_800Hz_noisy);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
x = data'; k = 2;
[m,v,w] = gaussmix2(x,[],[],k,'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 3;
[m2,v2,w2] = gaussmix2(x,[],[],k, 'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 4;
[m3,v3,w3] = gaussmix2(x,[],[],k, 'vf');    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% iteration, loop for the 10 different model orders
data = power_at_800Hz_noisy;
loop_AIC = 0; 
N = length(data);

%we find the AIC and MDL for each case using the following code
for q = 1 : 30 
    % We use GMM
    % Create an object with this GMM
    GMFIT = gmdistribution.fit(data', q);

    % We use the AIC
    loop_AIC(q,1) = GMFIT.AIC;
end

% plot the AIC with q, q = [1:10]
figure; set(gcf,'Color','w');
plot([1:30], loop_AIC);
grid on; grid minor;

% name the axes
set(gca,'FontSize',15);
title('AIC for GMM'); xlabel('Model Order p '); ylabel('AIC');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use k = 5
k = 5;
GMFIT = gmdistribution.fit(data', k);
%[m, v, w] = gaussmix(data', [], [], k);    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------












%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, all the active frames will be found using a VAD.

% We need to use the filename and change it every time
filename = 'clean_S_01_01';

% Initialize the variables for the for loop
store11 = [];

% Read the .wav file
[x, Srate, bits] = wavread( filename);	

% Initialize the matrix
total_noise_power_frame = [];

% Initialize the matrix
total_speech_active_frames_noisy = [];

% Define fs
Srate = 16000;

% Define the frame size in samples
% len = txinc_samples*10;
len = 320;

f=enframe(x,320,320/2);

% window overlap in percent of frame size
%PERC = 90;
PERC = 50;

len1 = floor(len * PERC/100);
len2 = len - len1;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the parameters
Nframes = floor(length(x)/len2)-1;
nFFT = 2 * len;

x_old = zeros(len1,1);
xfinal = zeros(Nframes*len2, 1);

noise_only_final = zeros(Nframes*len2, 1);
n_old = zeros(len1,1);

k = 1;
img = sqrt(-1);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

store_smoothed_p = [];
transient_presence_prob = 0.3;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define window
win = hamming(length(f(1,:)'));

%--- allocate memory and initialize various variables
k=1;
img=sqrt(-1);
x_old=zeros(len1,1);
Nframes=floor(length(x)/len2)-1;
xfinal=zeros(Nframes*len2,1);

% --------------- Initialize parameters ------------
k=1;
aa=0.98;
eta= 0.15;
mu=0.98;
c=sqrt(pi)/2;
qk=0.3;
qkr=(1-qk)/qk;
ksi_min=10^(-25/10); 

power_noise = [];
index_power = [];
index_power2 = [];

        % Main for loop
        for n = 1 : Nframes
            % Define the input signal
            insign = f(n,:)' .* win;

            % Take Fourier transform of  frame
            spec = fft(insign, nFFT);
        
            sig = abs(spec);
            sig2 = sig.^2;
            
            if (n == 1)
               noise_mu2 = sig2;
            end
            
            % posteriori SNR
            gammak=min(sig2./noise_mu2,40);  

            if n==1
                ksi=aa+(1-aa)*max(gammak-1,0);
            else
                ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);   

                % decision-direct estimate of a priori SNR
                ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
            end

            log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
            vad_decision = sum(log_sigma_k)/nFFT;    

            % Noise only frame found
            if (vad_decision < eta) 
                % store power 
                power_noise = [power_noise; spec'];

                % store the index
                index_power2 = [index_power2 n];
                
            else
                store11 = [store11; spec'];
                
            end

            vk=ksi.*gammak./(1+ksi);
            j0 = besseli(0,vk/2);
            j1 = besseli(1,vk/2);

            C=exp(-0.5*vk);
            A=((c*(vk.^0.5)).*C)./gammak;
            B=(1+vk).*j0+vk.*j1;
            hw=A.*B;

            sig=sig.*hw;

            % save for estimation of a priori SNR in next frame
            Xk_prev=sig.^2;  

            xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);

            xi_w= real( xi_w);

            xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
            x_old= xi_w(len1+ 1: len);

            k=k+len2;

            
        end
        
store11=store11(:,1:640/2+1);
store11 = real(store11);

% we use: store11
[afa fasfa] = size(store11);
variance_hat = [];
weight_hat = [];

for i = 1 : fasfa
    % the variance is the same in the real and imaginary parts    
    data = store11(:,i);

    data_mainData = data;
    clear data;

    nbins = 100;
    [y,x] = hist(data_mainData, nbins);
    [sigma, mu] = gaussfit_new(x, y);

%     data_mainData = data_mainData - mu;
%     [y,x] = hist(data_mainData, nbins);
%     [sigma, mu] = gaussfit_new(x, y);

    % in the end, we use "sigma" only
    clear mu;
    clear y;
    clear x;
    clear data;
    clear data_mainData;

    % we find variance 
    sigma_squared = sigma^2;

    % store the variance
    variance_hat = [variance_hat; sigma_squared];
end

% the variance is the same in the real and imaginary parts    
% we find the mean power parameter
variance_hat = 2 * variance_hat;

% the 3 parameters are: variance_hat, weight_hat

% we find the mean power parameter
new_mean_power_parameter2 = variance_hat';
new_mean_power_parameter2 = [new_mean_power_parameter2 flipud(new_mean_power_parameter2(2:end-1))];

% Store the variable
%save('new_clean_speech_stored_data2.mat', 'new_mean_power_parameter2');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------


















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, all the active frames will be found using a VAD.

% Clean speech signals will be utilised.
% The database "IEEE_corpus_clean_16kHz" will only be used.

% % We need to use the filename and change it every time
% filename = 'clean_S_01_01';
% 
% % Initialize the variables for the for loop
% total_speech_active_frames = [];
% 
% % Outter for loop
% for outter_loop = 1 : 10
%     % Inner for loop
%     for inner_loop = 1 : 72
%         % We change the filename based on outter_loop
%         % If statement for changing the filename based on outter_loop
%         if (outter_loop ~= 10)
%             filename(13) = num2str(outter_loop);
%         else
%             filename(12) = num2str(1);
%             filename(13) = num2str(0);
%         end
% 
%         % We change the filename based on k
%         % We use a temporary variable
%         temp_variable = num2str(inner_loop);
%     
%         % We now change the filename based on k
%         if (inner_loop < 10) 
%             filename(9) = num2str(0);
%             filename(10) = num2str(temp_variable(1));
%         else
%             filename(9) = num2str(temp_variable(1));
%             filename(10) = num2str(temp_variable(2));
%         end
% 
%         % Read the .wav file
%         [x, Srate, bits] = wavread( filename);	
% 
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%             % Speech and noise frame found
%             else
%                 total_speech_active_frames = [total_speech_active_frames; (x(k:k+len-1))'];
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
% 
%         end
%     end
% end
% 
% % Store the variable
% save('clean_speech_frames.mat', 'total_speech_active_frames');

% We use: "clean_speech_frames.mat".
voiced_clean_speech_frames = load('clean_speech_frames.mat');
% We use: "voiced_clean_speech_frames.total_speech_active_frames"

% We use: "size(voiced_clean_speech_frames.total_speech_active_frames)".
% The size is 177350 x 320.
% This means that 177350 frames with voiced speech exist. 

% Define the frame size in samples
len = size(voiced_clean_speech_frames.total_speech_active_frames,2); 

% Define the window
win = hanning(len)';  

% normalize window for equal level output 
win = win*len/sum(win);  

% Defien the matrix for the power at 800 Hz
power_at_800Hz = [];

% Define Fs
Fs = 16000;

% Define N
N = size(voiced_clean_speech_frames.total_speech_active_frames,2);

power_at_800Hz = [];

for i = 1 : size(voiced_clean_speech_frames.total_speech_active_frames,1)
    % We use the first frame: "voiced_clean_speech_frames.total_speech_active_frames(i,:)"
    insign = win .* voiced_clean_speech_frames.total_speech_active_frames(i,:);
    % We have used windowing.

    % Perform the FFT
    xdft = fft(insign, 2*length(win));
    
    power_at_800Hz = [power_at_800Hz; xdft];

%     % Perform the FFT
%     xdft = fft(insign);
%     xdft = xdft(1:N/2+1);
% 
%     % Find the power spectral density
%     psdx = (1/(Fs*N)) * abs(xdft).^2;
%     psdx(2:end-1) = 2*psdx(2:end-1);
% 
%     % Define the frequency axis
%     freq = 0:Fs/length(insign):Fs/2;
% 
%     % plot(freq,10*log10(psdx)); 
%     % grid on; title('Periodogram Using FFT');
%     % xlabel('Frequency (Hz)');
%     % ylabel('Power/Frequency (dB/Hz)');
% 
%     % Find the power at 800 Hz
%     idx = find(freq==800);
%     power_at_800Hz = [power_at_800Hz psdx(idx)];
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% plot the histogram
nbins = 10000;
[y,x] = hist(10*log10(power_at_800Hz), nbins);

% Plot the histogram
figure; set(gcf,'Color','w');
bar(x, y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the voiced clean speech power at 800 Hz. Log power is used.');
ylabel('Counts'); xlabel('Log power: 10 * log10(power)');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% The matrix "total_speech_active_frames.total_speech_active_frames" has all the speech active frames.
% We use the matrix "power_at_800Hz".

% Create a GMM with k mixtures
k = 2;

% Use the function "gaussmix" to fit a Gaussian mixture pdf to a set of data observations 
%     M(k,p)   Mixture means, one row per mixture. 
%     V(k,p)   Mixture variances, one row per mixture. 
%     W(k,1)   Mixture weights, one per mixture. The weights will sum to unity. 

% We use "gaussmix"
[m, v, w] = gaussmix((10*log10(x))', [], [], k);    

% Now, v is k x d
% We want v to be 1 x d x k
v_new = zeros(1, size(v,2), size(v,1));
v_new(1,:,1) = v(1,:);
v_new(1,:,2) = v(2,:);

% Create an object with this GMM
obj = gmdistribution(m, v_new, w);

% Plot the GMM
y = pdf(obj,(10*log10(x))');

%   LP(n,1) = log probability of each data point
%   RP(n,k) = relative probability of each mixture
%   KH(n,1) = highest probability mixture
%   KP(n,1) = relative probability of highest probability mixture
[LP, RP, KH, KP] = gaussmixp((10*log10(x))',m,v,w);


asfdda










%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% %iteration, loop for the 10 different model orders: ar(1) ... ar(10)
% loop_MDL = 0; 
% loop_AIC = 0; 
% Ep = 0; 
% N = length(10*log10(x));
% %coeffs = 0; 
% %estimate_matrix_vowel_a = 0;
% 
% %we find the AIC and MDL for each case using the following code
% for q = 1 : 10 
%     %loop_ar_coeffs = aryule(a_letter_samples, q);
%     
%     %we will use loop_ar_reflect_coeffs, which are the reflection coefficients
%     %we initialize the q first terms for the estimate of the data
%     %for w = 1 : q
%     %    estimate_matrix_vowel_a(w,1) = a_letter_samples(w,1);
%     %end
%     
%     %define coefficients
%     %coeffs = 0;
%     %for w = 1 : q
%     %    coeffs(w,1) = (-1) * loop_ar_coeffs(1,w+1);
%     %end
%     %we find the estimate of the data using ar(q)
%     %for i = (q+1) : N
%     %    %we find the estimate 
%     %    estimate_matrix_vowel_a(i,1) = dot(coeffs, estimate_matrix_vowel_a((1+i-(q+1)):(i-1),1));
%     %end
%     
%     
%     
%     %we find the loss function Ep, we use the "norm" command
%     Ep(q,1) = power(norm(a_letter_samples - estimate_matrix_vowel_a),2);
%     
%     %we use the MDL, where MDL = log(E) + (1/d(1,1)) * (p * log(d(1,1)));
%     loop_MDL(q,1) = log(Ep(q,1)) + (1/N) * (q * log(N));
%     
%     %we use the AIC, where AIC = log(E) + (2*p) / d(1,1);
%     loop_AIC(q,1) = log(Ep(q,1)) + (2*q) / N;
% end
% 
% %plot the MDL with q, q = [1:10]
% figure; set(gcf,'Color','w'); 
% subplot(1,2,1); plot([1:10], -loop_MDL);
% %name the axes
% set(gca,'FontSize',15);
% title('Letter a, N = 1000, MDL, AR(2) seems to be the optimal order'); xlabel('Model Order p '); ylabel('MDL');
% 
% %plot the AIC with q, q = [1:10]
% subplot(1,2,2); plot([1:10], -loop_AIC);
% %name the axes
% set(gca,'FontSize',15);
% title('AIC, We will also use AR(10) as over-modelling of data'); xlabel('Model Order p '); ylabel('AIC');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------














%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, all the active frames will be found using a VAD.

% Now, noisy speech signals will be utilised.
% The database "IEEE_corpus_sps_16kHz" will be used.

% % We need to use the filename and change it every time
% filename = 'ns_sps_0dB_S_01_01';
% 
% % Initialize the variables for the for loop
% total_speech_active_frames_noisy = [];
% 
% % Outter for loop
% for outter_loop = 1 : 10
%     % Inner for loop
%     for inner_loop = 1 : 72
%         % We change the filename based on outter_loop
%         % If statement for changing the filename based on outter_loop
%         if (outter_loop ~= 10)
%             filename(17) = num2str(0);
%             filename(18) = num2str(outter_loop);
%         else
%             filename(17) = num2str(1);
%             filename(18) = num2str(0);
%         end
% 
%         % We change the filename based on k
%         % We use a temporary variable
%         temp_variable = num2str(inner_loop);
%     
%         % We now change the filename based on k
%         if (inner_loop < 10) 
%             filename(14) = num2str(0);
%             filename(15) = num2str(temp_variable(1));
%         else
%             filename(14) = num2str(temp_variable(1));
%             filename(15) = num2str(temp_variable(2));
%         end
% 
%         % Read the .wav file
%         [x, Srate, bits] = wavread( filename);	
% 
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%             % Speech and noise frame found
%             else
%                 total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
% 
%         end
%     end
% end
%
% % Store the variable
% save('noisy_1_speech_frames.mat', 'total_speech_active_frames_noisy');

% We use: "noisy_1_speech_frames.mat".
noisy_1_speech_frames = load('noisy_1_speech_frames');
% We use: noisy_1_speech_frames.total_speech_active_frames_noisy
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % We now use the "Noisy_Speech_Database_white_noise_16kHz" database.
% filename = 'sp01_white_noise_0dB';
% %filename = 'sp01_white_noise_5dB';
% 
% % Initialize the variables for the for loop
% total_speech_active_frames_noisy = [];
% 
% % Outter for loop
% for outter_loop = 1 : 30
%     for inner_loop = 1 : 2
%         % We change the filename based on outter_loop
%         % We use a temporary variable
%         temp_variable = num2str(outter_loop);
% 
%         % We now change the filename based on k
%         if (outter_loop < 10) 
%             filename(3) = num2str(0);
%             filename(4) = num2str(temp_variable(1));
%         else
%             filename(3) = num2str(temp_variable(1));
%             filename(4) = num2str(temp_variable(2));
%         end
% 
%         % If statement to change the filename
%         if (inner_loop == 1) 
%             filename(18) = num2str(0);
%         else 
%             filename(18) = num2str(5);
%         end
%         
%         if (outter_loop ~= 21 && outter_loop ~= 24)
%             % Read the .wav file
%             [x, Srate, bits] = wavread( filename);	
% 
%             % Define the frame size in samples
%             len=floor(20*Srate/1000); 
%             if rem(len,2)==1, len=len+1; end;
% 
%             % window overlap in percent of frame size
%             PERC=50; 
%             len1=floor(len*PERC/100);
%             len2=len-len1;
% 
%             % define window
%             win=hanning(len);  
% 
%             % normalize window for equal level output 
%             win = win*len2/sum(win);  
% 
%             % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%             nFFT=2*len;
%             j=1;
%             noise_mean=zeros(nFFT,1);
%             for k=1:6
%                 noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%                 j=j+len;
%             end
%             noise_mu=noise_mean/6;
%             noise_mu2=noise_mu.^2;
% 
%             %--- allocate memory and initialize various variables
%             k=1;
%             img=sqrt(-1);
%             x_old=zeros(len1,1);
%             Nframes=floor(length(x)/len2)-1;
%             xfinal=zeros(Nframes*len2,1);
% 
%             % --------------- Initialize parameters ------------
%             k=1;
%             aa=0.98;
%             eta= 0.15;
%             mu=0.98;
%             c=sqrt(pi)/2;
%             qk=0.3;
%             qkr=(1-qk)/qk;
%             ksi_min=10^(-25/10); 
% 
%             % Main for loop
%             for n=1:Nframes
% 
%                 insign=win.*x(k:k+len-1);
% 
%                 %--- Take fourier transform of  frame
%                 %
%                 spec=fft(insign,nFFT);
%                 sig=abs(spec); % compute the magnitude
%                 sig2=sig.^2;
% 
%                 gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%                 if n==1
%                     ksi=aa+(1-aa)*max(gammak-1,0);
%                 else
%                     ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                     % decision-direct estimate of a priori SNR
%                     ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%                 end
% 
%                 log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%                 vad_decision = sum(log_sigma_k)/nFFT;    
% 
%                 % Noise only frame found
%                 if (vad_decision < eta) 
%                     noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%                 % Speech and noise frame found
%                 else
%                     total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
%                 end
% 
%                 vk=ksi.*gammak./(1+ksi);
%                 j0 = besseli(0,vk/2);
%                 j1 = besseli(1,vk/2);
% 
%                 C=exp(-0.5*vk);
%                 A=((c*(vk.^0.5)).*C)./gammak;
%                 B=(1+vk).*j0+vk.*j1;
%                 hw=A.*B;
% 
%                 sig=sig.*hw;
% 
%                 % save for estimation of a priori SNR in next frame
%                 Xk_prev=sig.^2;  
% 
%                 xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%                 xi_w= real( xi_w);
% 
%                 xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%                 x_old= xi_w(len1+ 1: len);
% 
%                 k=k+len2;
% 
%             end
%         end 
%     end
% end
% 
% filename = 'sp01_white_noise_10dB';
% %filename = 'sp01_white_noise_15dB';
% 
% % Outter for loop
% for outter_loop = 1 : 30
%     for inner_loop = 1 : 2
%         % We change the filename based on i
%         % We use a temporary variable
%         temp_variable = num2str(outter_loop);
% 
%         % We now change the filename based on k
%         if (outter_loop < 10) 
%             filename(3) = num2str(0);
%             filename(4) = num2str(temp_variable(1));
%         else
%             filename(3) = num2str(temp_variable(1));
%             filename(4) = num2str(temp_variable(2));
%         end
% 
%         % If statement to change the filename
%         if (inner_loop == 1) 
%             filename(19) = num2str(0);
%         else 
%             filename(19) = num2str(5);
%         end
%         
%         if (outter_loop ~= 21 && outter_loop ~= 24)
%             % Read the .wav file
%             [x, Srate, bits] = wavread( filename);	
% 
%             % Define the frame size in samples
%             len=floor(20*Srate/1000); 
%             if rem(len,2)==1, len=len+1; end;
% 
%             % window overlap in percent of frame size
%             PERC=50; 
%             len1=floor(len*PERC/100);
%             len2=len-len1;
% 
%             % define window
%             win=hanning(len);  
% 
%             % normalize window for equal level output 
%             win = win*len2/sum(win);  
% 
%             % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%             nFFT=2*len;
%             j=1;
%             noise_mean=zeros(nFFT,1);
%             for k=1:6
%                 noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%                 j=j+len;
%             end
%             noise_mu=noise_mean/6;
%             noise_mu2=noise_mu.^2;
% 
%             %--- allocate memory and initialize various variables
%             k=1;
%             img=sqrt(-1);
%             x_old=zeros(len1,1);
%             Nframes=floor(length(x)/len2)-1;
%             xfinal=zeros(Nframes*len2,1);
% 
%             % --------------- Initialize parameters ------------
%             k=1;
%             aa=0.98;
%             eta= 0.15;
%             mu=0.98;
%             c=sqrt(pi)/2;
%             qk=0.3;
%             qkr=(1-qk)/qk;
%             ksi_min=10^(-25/10); 
% 
%             % Main for loop
%             for n=1:Nframes
% 
%                 insign=win.*x(k:k+len-1);
% 
%                 %--- Take fourier transform of  frame
%                 %
%                 spec=fft(insign,nFFT);
%                 sig=abs(spec); % compute the magnitude
%                 sig2=sig.^2;
% 
%                 gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%                 if n==1
%                     ksi=aa+(1-aa)*max(gammak-1,0);
%                 else
%                     ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                     % decision-direct estimate of a priori SNR
%                     ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%                 end
% 
%                 log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%                 vad_decision = sum(log_sigma_k)/nFFT;    
% 
%                 % Noise only frame found
%                 if (vad_decision < eta) 
%                     noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%                 % Speech and noise frame found
%                 else
%                     total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
%                 end
% 
%                 vk=ksi.*gammak./(1+ksi);
%                 j0 = besseli(0,vk/2);
%                 j1 = besseli(1,vk/2);
% 
%                 C=exp(-0.5*vk);
%                 A=((c*(vk.^0.5)).*C)./gammak;
%                 B=(1+vk).*j0+vk.*j1;
%                 hw=A.*B;
% 
%                 sig=sig.*hw;
% 
%                 % save for estimation of a priori SNR in next frame
%                 Xk_prev=sig.^2;  
% 
%                 xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%                 xi_w= real( xi_w);
% 
%                 xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%                 x_old= xi_w(len1+ 1: len);
% 
%                 k=k+len2;
% 
%             end
%         end
%     end
% end
%
% % Store the variable
% save('noisy_2_speech_frames.mat', 'total_speech_active_frames_noisy');

% We use: "noisy_2_speech_frames.mat".
noisy_2_speech_frames = load('noisy_2_speech_frames');
% We use: noisy_2_speech_frames.total_speech_active_frames_noisy
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We now use the "NOIZEUS_Noisy_Speech_Signal_Database_16kHz" database.

% % Define the file name
% filename = 'sp01_train_sn0';
% 
% % We initialize the matrix.
% total_speech_active_frames_noisy = [];
% 
% % Outter for loop
% for outter_loop = 1 : 4
%     for inner_loop = 1 : 30
%         % We change the filename based on inner_loop
%         % We use a temporary variable
%         temp_variable = num2str(inner_loop);
% 
%         % We now change the filename based on inner_loop
%         if (inner_loop < 10) 
%             filename(3) = num2str(0);
%             filename(4) = num2str(temp_variable(1));
%         else
%             filename(3) = num2str(temp_variable(1));
%             filename(4) = num2str(temp_variable(2));
%         end
% 
%         % If statement to change the filename
%         if (outter_loop == 1) 
%             filename(14) = num2str(0);
%         elseif (outter_loop == 2) 
%             filename(14) = num2str(5);
%         elseif (outter_loop == 3) 
%             filename(14) = num2str(1);
%             filename(15) = num2str(0);
%         else 
%             filename(14) = num2str(1);
%             filename(15) = num2str(5);
%         end
% 
%         %if (strcmp(filename, 'sp04_babble_sn10') ~= 1)
%         % Read the .wav file
%         [x, Srate, bits] = wavread( filename);	
% 
%         % Define the frame size in samples
%         len=floor(20*Srate/1000); 
%         if rem(len,2)==1, len=len+1; end;
% 
%         % window overlap in percent of frame size
%         PERC=50; 
%         len1=floor(len*PERC/100);
%         len2=len-len1;
% 
%         % define window
%         win=hanning(len);  
% 
%         % normalize window for equal level output 
%         win = win*len2/sum(win);  
% 
%         % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
%         nFFT=2*len;
%         j=1;
%         noise_mean=zeros(nFFT,1);
%         for k=1:6
%             noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
%             j=j+len;
%         end
%         noise_mu=noise_mean/6;
%         noise_mu2=noise_mu.^2;
% 
%         %--- allocate memory and initialize various variables
%         k=1;
%         img=sqrt(-1);
%         x_old=zeros(len1,1);
%         Nframes=floor(length(x)/len2)-1;
%         xfinal=zeros(Nframes*len2,1);
% 
%         % --------------- Initialize parameters ------------
%         k=1;
%         aa=0.98;
%         eta= 0.15;
%         mu=0.98;
%         c=sqrt(pi)/2;
%         qk=0.3;
%         qkr=(1-qk)/qk;
%         ksi_min=10^(-25/10); 
% 
%         % Main for loop
%         for n=1:Nframes
% 
%             insign=win.*x(k:k+len-1);
% 
%             %--- Take fourier transform of  frame
%             %
%             spec=fft(insign,nFFT);
%             sig=abs(spec); % compute the magnitude
%             sig2=sig.^2;
% 
%             gammak=min(sig2./noise_mu2,40);  % posteriori SNR
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
%             % Speech and noise frame found
%             else
%                 total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
%         %end
%         end
%     end
% end
% 
% % Store the variable
% save('noisy_10_speech_frames.mat', 'total_speech_active_frames_noisy');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Load the frames

% We use: "noisy_3_speech_frames.mat".
noisy_3_speech_frames = load('noisy_3_speech_frames');
% We use: noisy_3_speech_frames.total_speech_active_frames_noisy

% We use: "noisy_4_speech_frames.mat".
noisy_4_speech_frames = load('noisy_4_speech_frames');
% We use: noisy_4_speech_frames.total_speech_active_frames_noisy

% We use: "noisy_5_speech_frames.mat".
noisy_5_speech_frames = load('noisy_5_speech_frames');
% We use: noisy_5_speech_frames.total_speech_active_frames_noisy

% We use: "noisy_6_speech_frames.mat".
noisy_6_speech_frames = load('noisy_6_speech_frames');
% We use: noisy_6_speech_frames.total_speech_active_frames_noisy

% We use: "noisy_7_speech_frames.mat".
noisy_7_speech_frames = load('noisy_7_speech_frames');
% We use: noisy_7_speech_frames.total_speech_active_frames_noisy

% We use: "noisy_8_speech_frames.mat".
noisy_8_speech_frames = load('noisy_8_speech_frames');
% We use: noisy_8_speech_frames.total_speech_active_frames_noisy

% We use: "noisy_9_speech_frames.mat".
noisy_9_speech_frames = load('noisy_9_speech_frames');
% We use: noisy_9_speech_frames.total_speech_active_frames_noisy

% We use: "noisy_10_speech_frames.mat".
noisy_10_speech_frames = load('noisy_10_speech_frames');
% We use: noisy_10_speech_frames.total_speech_active_frames_noisy
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Initialize the matrix with the noisy speech frames
final_noisy_frames = [];

% Define the matrix with the noisy speech frames
final_noisy_frames = [noisy_1_speech_frames.total_speech_active_frames_noisy; ...
    noisy_2_speech_frames.total_speech_active_frames_noisy; ...
    noisy_3_speech_frames.total_speech_active_frames_noisy; ...
    noisy_4_speech_frames.total_speech_active_frames_noisy; ...
    noisy_5_speech_frames.total_speech_active_frames_noisy; ...
    noisy_6_speech_frames.total_speech_active_frames_noisy; ...
    noisy_7_speech_frames.total_speech_active_frames_noisy; ...
    noisy_8_speech_frames.total_speech_active_frames_noisy; ...
    noisy_9_speech_frames.total_speech_active_frames_noisy; ...
    noisy_10_speech_frames.total_speech_active_frames_noisy];

% Define the frame size in samples
len = size(final_noisy_frames,2); 

% Define the window
win = hanning(len)';  

% normalize window for equal level output 
win = win*len/sum(win);  

% Define the matrix for the power at 800 Hz
power_at_800Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define N
N = size(final_noisy_frames,2);

% for loop for all frames
for i = 1 : size(final_noisy_frames,1)
    % We use the first frame: "voiced_clean_speech_frames.total_speech_active_frames(i,:)"
    insign = win .* final_noisy_frames(i,:);
    % We have used windowing.

    % Perform the FFT
    xdft = fft(insign);
    xdft = xdft(1:N/2+1);

    % Find the power spectral density
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    % Define the frequency axis
    freq = 0:Fs/length(insign):Fs/2;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the speech level in units of power
    lev = activlev(final_noisy_frames(i,:), Fs);     
    
    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)/lev];
end

% plot the histogram
nbins = 1000;
[y,x] = hist(log10(power_at_800Hz_noisy), nbins);

% Plot the histogram
figure; set(gcf,'Color','w');
bar(x, y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the voiced noisy speech power at 800 Hz. Log power is used.');
ylabel('Counts'); xlabel('Log power: 10 * log10(power)');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Fit a distribution to the data
%data = 10*log10(x);
data = 10*log10(power_at_800Hz_noisy);

% Compute and plot results
[D PD] = allfitdist(data,'PDF'); 

% Show output from best fit
% We use D(1)
D(1)
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the variable
save('gaussianPDF.mat', 'data');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Create a GMM model for the data

% Load data
data_GMM = load('gaussianPDF.mat');
data = data_GMM.data;

% create GMM with k mixtures
x = data'; k = 2;
[m,v,w] = gaussmix(x,[],[],k);    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 3;
[m2,v2,w2] = gaussmix(x,[],[],k);    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create GMM with k mixtures
k = 4;
[m3,v3,w3] = gaussmix(x,[],[],k);    
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------














%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%for loop = 1 : 10
sp = total_snr15{1}(:,1);
Fs = 16000;
fx = fxpefac(sp, Fs);

% We initialize the matrix.
total_speech_active_frames_noisy = [];

% Define a frame counter
frame_counter = [];

% Define the parameters
x = sp;
Srate = Fs;

% Define the frame size in samples
len=floor(20*Srate/1000); 
if rem(len,2)==1, len=len+1; end;

% window overlap in percent of frame size
PERC=50; 
len1=floor(len*PERC/100);
len2=len-len1;

% define window
win=hanning(len);  

% normalize window for equal level output 
win = win*len2/sum(win);  

% Noise magnitude calculations - assuming that the first 6 frames is noise/silence
nFFT=2*len;
j=1;
noise_mean=zeros(nFFT,1);
for k=1:6
    noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
    j=j+len;
end
noise_mu=noise_mean/6;
noise_mu2=noise_mu.^2;

%--- allocate memory and initialize various variables
k=1;
img=sqrt(-1);
x_old=zeros(len1,1);
Nframes=floor(length(x)/len2)-1;
xfinal=zeros(Nframes*len2,1);

% --------------- Initialize parameters ------------
k=1;
aa=0.98;
eta= 0.15;
mu=0.98;
c=sqrt(pi)/2;
qk=0.3;
qkr=(1-qk)/qk;
ksi_min=10^(-25/10); 

% Main for loop
for n=1:Nframes

    insign=win.*x(k:k+len-1);

    %--- Take fourier transform of  frame
    %
    spec=fft(insign,nFFT);
    sig=abs(spec); % compute the magnitude
    sig2=sig.^2;

    gammak=min(sig2./noise_mu2,40);  % posteriori SNR
    if n==1
        ksi=aa+(1-aa)*max(gammak-1,0);
    else
        ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
        % decision-direct estimate of a priori SNR
        ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
    end

    log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
    vad_decision = sum(log_sigma_k)/nFFT;    

    % Noise only frame found
    if (vad_decision < eta) 
        noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
    % Speech and noise frame found
    else
        total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
        frame_counter = [frame_counter n];
    end

    vk=ksi.*gammak./(1+ksi);
    j0 = besseli(0,vk/2);
    j1 = besseli(1,vk/2);

    C=exp(-0.5*vk);
    A=((c*(vk.^0.5)).*C)./gammak;
    B=(1+vk).*j0+vk.*j1;
    hw=A.*B;

    sig=sig.*hw;

    % save for estimation of a priori SNR in next frame
    Xk_prev=sig.^2;  

    xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);

    xi_w= real( xi_w);

    xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
    x_old= xi_w(len1+ 1: len);

    k=k+len2;
end


% We use the frame counter
frame_counter = frame_counter(frame_counter <= length(fx));
fx = fx(frame_counter);



% Define the matrix with the noisy speech frames
final_noisy_frames = total_speech_active_frames_noisy;

% Define the frame size in samples
len = size(final_noisy_frames,2); 

% Define the window
win = hanning(len)';  

% normalize window for equal level output 
win = win*len/sum(win);  

% Define the matrix for the power at 800 Hz
power_at_800Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define N
N = size(final_noisy_frames,2);

for i = 1 : size(final_noisy_frames,1)
    % We use the first frame: "voiced_clean_speech_frames.total_speech_active_frames(i,:)"
    insign = win .* final_noisy_frames(i,:);
    % We have used windowing.

    % Perform the FFT
    xdft = fft(insign);
    xdft = xdft(1:N/2+1);

    % Find the power spectral density
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    psdx = psdx ./ (fx(i)^2);
    
    % Define the frequency axis
    freq = 0:Fs/length(insign):Fs/2;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)];
end

% plot the histogram
nbins = 100000;
[y,x] = hist(power_at_800Hz_noisy, nbins);

% Plot the histogram
figure; set(gcf,'Color','w');
bar(10*log10(x), y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the voiced noisy speech power at 800 Hz. Log power is used.');
ylabel('Counts'); xlabel('Log power: 10 * log10(power)');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------














%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % Define the filename
% filename = 'sp01_airport_sn15';
% 
% % Read the .wav file
% [sp, Fs, bits] = wavread( filename);	


% sp = total_snr15{1}(:,1);
% Fs = 16000;

% We now use the "NOIZEUS_Noisy_Speech_Signal_Database_16kHz" database.
%filename = 'sp01_airport_sn0';
filename = 'sp01_airport_sn5';
%filename = 'sp01_airport_sn10';
%filename = 'sp01_airport_sn15';

% We initialize the matrix.
total_speech_active_frames_noisy = [];

% Outter for loop
for outter_loop = 1 : 4
    for inner_loop = 1 : 30
        % We change the filename based on inner_loop
        % We use a temporary variable
        temp_variable = num2str(inner_loop);

        % We now change the filename based on inner_loop
        if (inner_loop < 10) 
            filename(3) = num2str(0);
            filename(4) = num2str(temp_variable(1));
        else
            filename(3) = num2str(temp_variable(1));
            filename(4) = num2str(temp_variable(2));
        end

        % If statement to change the filename
        if (outter_loop == 1) 
            filename(16) = num2str(0);
        elseif (outter_loop == 2) 
            filename(16) = num2str(5);
        elseif (outter_loop == 3) 
            filename(16) = num2str(1);
            filename(17) = num2str(0);
        else 
            filename(16) = num2str(1);
            filename(17) = num2str(5);
        end

        % Use the MMSE algorithm
        outfile = strcat(filename, '_out');
        mmse(filename, outfile, 0);
        
        % Read the .wav file
        [x, Srate, bits] = wavread(outfile);	

        % Define the frame size in samples
        len=floor(20*Srate/1000); 
        if rem(len,2)==1, len=len+1; end;

        % window overlap in percent of frame size
        PERC=50; 
        len1=floor(len*PERC/100);
        len2=len-len1;

        % define window
        win=hanning(len);  

        % normalize window for equal level output 
        win = win*len2/sum(win);  

        % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
        nFFT=2*len;
        j=1;
        noise_mean=zeros(nFFT,1);
        for k=1:6
            noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
            j=j+len;
        end
        noise_mu=noise_mean/6;
        noise_mu2=noise_mu.^2;

        %--- allocate memory and initialize various variables
        k=1;
        img=sqrt(-1);
        x_old=zeros(len1,1);
        Nframes=floor(length(x)/len2)-1;
        xfinal=zeros(Nframes*len2,1);

        % --------------- Initialize parameters ------------
        k=1;
        aa=0.98;
        eta= 0.15;
        mu=0.98;
        c=sqrt(pi)/2;
        qk=0.3;
        qkr=(1-qk)/qk;
        ksi_min=10^(-25/10); 

        % Main for loop
        for n=1:Nframes

            insign=win.*x(k:k+len-1);

            %--- Take fourier transform of  frame
            %
            spec=fft(insign,nFFT);
            sig=abs(spec); % compute the magnitude
            sig2=sig.^2;

            gammak=min(sig2./noise_mu2,40);  % posteriori SNR
            if n==1
                ksi=aa+(1-aa)*max(gammak-1,0);
            else
                ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
                % decision-direct estimate of a priori SNR
                ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
            end

            log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
            vad_decision = sum(log_sigma_k)/nFFT;    

            % Noise only frame found
            if (vad_decision < eta) 
                noise_mu2= mu* noise_mu2+ (1- mu)* sig2;
            % Speech and noise frame found
            else
                total_speech_active_frames_noisy = [total_speech_active_frames_noisy; (x(k:k+len-1))'];
            end

            vk=ksi.*gammak./(1+ksi);
            j0 = besseli(0,vk/2);
            j1 = besseli(1,vk/2);

            C=exp(-0.5*vk);
            A=((c*(vk.^0.5)).*C)./gammak;
            B=(1+vk).*j0+vk.*j1;
            hw=A.*B;

            sig=sig.*hw;

            % save for estimation of a priori SNR in next frame
            Xk_prev=sig.^2;  

            xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);

            xi_w= real( xi_w);

            xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
            x_old= xi_w(len1+ 1: len);

            k=k+len2;

        end
       
    end
end

% Store the variable
save('afterMMSE_voice_speech_signal_2.mat', 'total_speech_active_frames_noisy');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------





total_speech_active_frames_noisy2 = load('afterMMSE_voice_speech_signal_2.mat');
total_speech_active_frames_noisy3 = load('afterMMSE_voice_speech_signal_3.mat');


final_noisy_frames = [total_speech_active_frames_noisy2.total_speech_active_frames_noisy; total_speech_active_frames_noisy3.total_speech_active_frames_noisy];

% Define the frame size in samples
len = size(final_noisy_frames,2); 

% Define the window
win = hanning(len)';  

% normalize window for equal level output 
win = win*len/sum(win);  

% Define the matrix for the power at 800 Hz
power_at_800Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define N
N = size(final_noisy_frames,2);

for i = 1 : size(final_noisy_frames,1)
    % We use the first frame: "voiced_clean_speech_frames.total_speech_active_frames(i,:)"
    insign = win .* final_noisy_frames(i,:);
    % We have used windowing.

    % Perform the FFT
    xdft = fft(insign);
    xdft = xdft(1:N/2+1);

    % Find the power spectral density
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    % Define the frequency axis
    freq = 0:Fs/length(insign):Fs/2;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)];
end




power_at_800Hz_noisy2 = power_at_800Hz_noisy;

% Initialize the matrix with the noisy speech frames
final_noisy_frames = [noisy_2_speech_frames.total_speech_active_frames_noisy; noisy_3_speech_frames.total_speech_active_frames_noisy];

% Define the frame size in samples
len = size(final_noisy_frames,2); 

% Define the window
win = hanning(len)';  

% normalize window for equal level output 
win = win*len/sum(win);  

% Define the matrix for the power at 800 Hz
power_at_800Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define N
N = size(final_noisy_frames,2);

for i = 1 : size(final_noisy_frames,1)
    % We use the first frame: "voiced_clean_speech_frames.total_speech_active_frames(i,:)"
    insign = win .* final_noisy_frames(i,:);
    % We have used windowing.

    % Perform the FFT
    xdft = fft(insign);
    xdft = xdft(1:N/2+1);

    % Find the power spectral density
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    % Define the frequency axis
    freq = 0:Fs/length(insign):Fs/2;

    % plot(freq,10*log10(psdx)); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy psdx(idx)];
end

% plot the histogram
nbins = 100000;
[y,x] = hist(power_at_800Hz_noisy/power_at_800Hz_noisy2(1:length(power_at_800Hz_noisy)), nbins);

% Plot the histogram
figure; set(gcf,'Color','w');
bar(10*log10(x), y, 'hist');

% Set the title of the image
set(gca,'FontSize',15);
title('Histogram for the voiced noisy speech power at 800 Hz. Log power is used.');
ylabel('Counts'); xlabel('Log power: 10 * log10(power)');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------


















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use the matrix "power_at_800Hz".

% Create a GMM with k mixtures
k = 2;

% Use the function "gaussmix" to fit a Gaussian mixture pdf to a set of data observations 
%     M(k,p)   Mixture means, one row per mixture. 
%     V(k,p)   Mixture variances, one row per mixture. 
%     W(k,1)   Mixture weights, one per mixture. The weights will sum to unity. 

% We use the "gaussmix"
[m, v, w] = gaussmix((10*log10(x))', [], [], k);    

% Now, v is k x d
% We want v to be 1 x d x k
v_new = zeros(1, size(v,2), size(v,1));
v_new(1,:,1) = v(1,:);
v_new(1,:,2) = v(2,:);

% Create an object with this GMM
obj = gmdistribution(m, v_new, w);

%   LP(n,1) = log probability of each data point
%   RP(n,k) = relative probability of each mixture
%   KH(n,1) = highest probability mixture
%   KP(n,1) = relative probability of highest probability mixture
[LP, RP, KH, KP] = gaussmixp(total_speech_active_frames,m,v,w);













%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this section of the code, we apply the EM algorithm for the GMM.
% We apply the GMM algorithm: we use the "emgmm.m" from the "stprtool" folder.

% define the path for the toolbox
stprpath;

% We define the "total array" variable.
% We use the previous matrix for the total speech active frames.
total_array = 10*log10(power_at_800Hz_noisy);

% We use K = 4 Gaussian components (ncomp)
% We display the information (we use: verb)
model = emgmm(double(total_array), struct('ncomp', 4, 'verb', 1));
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% %the model has: "Cov: [AxAx4 double]"
% 
% %display the Covariance images
% figure; set(gcf,'Color','w');
% %display the first Covariance image
% subplot(1,4,1); 
% imagesc(model.Cov(:,:,1));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 5: 1000 Column Vectors. Covariance Matrices.');
% 
% %display the second Covariance image
% subplot(1,4,2); 
% imagesc(model.Cov(:,:,2));
% 
% %display the third Covariance image
% subplot(1,4,3); 
% imagesc(model.Cov(:,:,3));
% 
% %display the fourth Covariance image
% subplot(1,4,4); 
% imagesc(model.Cov(:,:,4));
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% %the model has: "Mean: [Ax4 double]"
% 
% %display the mean images
% figure; set(gcf,'Color','w');
% %display the first mean image
% subplot(1,4,1); 
% imagesc(imresize(model.Mean(:,1), [5 5], 'bilinear'));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 5: 1000 Column Vectors. Mean Matrices.');
% 
% %display the second mean image
% subplot(1,4,2); 
% imagesc(imresize(model.Mean(:,2), [5 5], 'bilinear'));
% 
% %display the third mean image
% subplot(1,4,3); 
% imagesc(imresize(model.Mean(:,3), [5 5], 'bilinear'));
% 
% %display the fourth mean image
% subplot(1,4,4); 
% imagesc(imresize(model.Mean(:,4), [5 5], 'bilinear'));
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We re-apply the GMM algorithm: we use the "emgmm.m".
% We now use random initialisation and full covariance.
model_random_full = emgmm(double(total_array), struct('ncomp',4, 'verb',1, 'init','random', 'cov_type','full'));

% We now use cmeans initialisation and full covariance.
model_cmeans_full = emgmm(double(total_array), struct('ncomp',4, 'verb',1, 'init','cmeans', 'cov_type','full'));

% We now use random initialisation and diagonal covariance.
model_random_diag = emgmm(double(total_array), struct('ncomp',4, 'verb',1, 'init','random', 'cov_type','diag'));

% We now use cmeans initialisation and diagonal covariance.
model_cmeans_diag = emgmm(double(total_array), struct('ncomp',4, 'verb',1, 'init','cmeans', 'cov_type','diag'));
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% %we use "model_random_full"
% 
% %display the mean images
% figure; set(gcf,'Color','w');
% %display the first mean image
% subplot(2,4,1); 
% imagesc(imresize(model_random_full.Mean(:,1), [5 5], 'bilinear'));
% 
% %display the second mean image
% subplot(2,4,2); 
% imagesc(imresize(model_random_full.Mean(:,2), [5 5], 'bilinear'));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Random Initialisation and Full Covariance. Mean Matrices.');
% 
% %display the third mean image
% subplot(2,4,3); 
% imagesc(imresize(model_random_full.Mean(:,3), [5 5], 'bilinear'));
% 
% %display the fourth mean image
% subplot(2,4,4); 
% imagesc(imresize(model_random_full.Mean(:,4), [5 5], 'bilinear'));
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% %we use "model_random_full"
% %display the Covariance images
% 
% %display the first Covariance image
% subplot(2,4,5); 
% imagesc(model_random_full.Cov(:,:,1));
% 
% %display the second Covariance image
% subplot(2,4,6); 
% imagesc(model_random_full.Cov(:,:,2));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Random Initialisation and Full Covariance. Covariance Matrices.');
% 
% %display the third Covariance image
% subplot(2,4,7); 
% imagesc(model_random_full.Cov(:,:,3));
% 
% %display the fourth Covariance image
% subplot(2,4,8); 
% imagesc(model_random_full.Cov(:,:,4));
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% %we use "model_cmeans_full"
% 
% %display the mean images
% figure; set(gcf,'Color','w');
% %display the first mean image
% subplot(2,4,1); 
% imagesc(imresize(model_cmeans_full.Mean(:,1), [5 5], 'bilinear'));
% 
% %display the second mean image
% subplot(2,4,2); 
% imagesc(imresize(model_cmeans_full.Mean(:,2), [5 5], 'bilinear'));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Cmeans Initialisation and Full Covariance. Mean Matrices.');
% 
% %display the third mean image
% subplot(2,4,3); 
% imagesc(imresize(model_cmeans_full.Mean(:,3), [5 5], 'bilinear'));
% 
% %display the fourth mean image
% subplot(2,4,4); 
% imagesc(imresize(model_cmeans_full.Mean(:,4), [5 5], 'bilinear'));
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% %we use "model_cmeans_full"
% %display the Covariance images
% 
% %display the first Covariance image
% subplot(2,4,5); 
% imagesc(model_cmeans_full.Cov(:,:,1));
% 
% %display the second Covariance image
% subplot(2,4,6); 
% imagesc(model_cmeans_full.Cov(:,:,2));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Cmeans Initialisation and Full Covariance. Covariance Matrices.');
% 
% %display the third Covariance image
% subplot(2,4,7); 
% imagesc(model_cmeans_full.Cov(:,:,3));
% 
% %display the fourth Covariance image
% subplot(2,4,8); 
% imagesc(model_cmeans_full.Cov(:,:,4));
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% %we use "model_random_diag"
% 
% %display the mean images
% figure; set(gcf,'Color','w');
% %display the first mean image
% subplot(2,4,1); 
% imagesc(imresize(model_random_diag.Mean(:,1), [5 5], 'bilinear'));
% 
% %display the second mean image
% subplot(2,4,2); 
% imagesc(imresize(model_random_diag.Mean(:,2), [5 5], 'bilinear'));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Random Initialisation and Diagonal Covariance. Mean Matrices.');
% 
% %display the third mean image
% subplot(2,4,3); 
% imagesc(imresize(model_random_diag.Mean(:,3), [5 5], 'bilinear'));
% 
% %display the fourth mean image
% subplot(2,4,4); 
% imagesc(imresize(model_random_diag.Mean(:,4), [5 5], 'bilinear'));
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% %we use "model_random_diag"
% %display the Covariance images
% 
% %display the first Covariance image
% subplot(2,4,5); 
% imagesc(model_random_diag.Cov(:,:,1));
% 
% %display the second Covariance image
% subplot(2,4,6); 
% imagesc(model_random_diag.Cov(:,:,2));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Random Initialisation and Diagonal Covariance. Covariance Matrices.');
% 
% %display the third Covariance image
% subplot(2,4,7); 
% imagesc(model_random_diag.Cov(:,:,3));
% 
% %display the fourth Covariance image
% subplot(2,4,8); 
% imagesc(model_random_diag.Cov(:,:,4));
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% %we use "model_cmeans_diag"
% 
% %display the mean images
% figure; set(gcf,'Color','w');
% %display the first mean image
% subplot(2,4,1); 
% imagesc(imresize(model_cmeans_diag.Mean(:,1), [5 5], 'bilinear'));
% 
% %display the second mean image
% subplot(2,4,2); 
% imagesc(imresize(model_cmeans_diag.Mean(:,2), [5 5], 'bilinear'));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Cmeans Initialisation and Diagonal Covariance. Mean Matrices.');
% 
% %display the third mean image
% subplot(2,4,3); 
% imagesc(imresize(model_cmeans_diag.Mean(:,3), [5 5], 'bilinear'));
% 
% %display the fourth mean image
% subplot(2,4,4); 
% imagesc(imresize(model_cmeans_diag.Mean(:,4), [5 5], 'bilinear'));
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% %we use "model_cmeans_diag"
% %display the Covariance images
% 
% %display the first Covariance image
% subplot(2,4,5); 
% imagesc(model_cmeans_diag.Cov(:,:,1));
% 
% %display the second Covariance image
% subplot(2,4,6); 
% imagesc(model_cmeans_diag.Cov(:,:,2));
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Cmeans Initialisation and Diagonal Covariance. Covariance Matrices.');
% 
% %display the third Covariance image
% subplot(2,4,7); 
% imagesc(model_cmeans_diag.Cov(:,:,3));
% 
% %display the fourth Covariance image
% subplot(2,4,8); 
% imagesc(model_cmeans_diag.Cov(:,:,4));
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% %plot the Log-Likelihoods
% figure; set(gcf,'Color','w');
% plot(model_random_full.logL, 'b'); hold on;
% plot(model_cmeans_full.logL, 'm'); hold on;
% plot(model_random_diag.logL, 'g'); hold on;
% plot(model_cmeans_diag.logL, 'r'); hold off;
% 
% 
% %set the title of the image
% set(gca,'FontSize',15);
% title('Exercise 6: Comparison of Log-Likelihoods.');
% %set the labels and the legend
% xlabel('Iterations'); ylabel('Log - Likelihood');
% legend('Random and Full', 'Cmeans and Full', 'Random and Diag', 'Cmeans and Diag');
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------


















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, the voicebox will be used. 
% Specifically, the "fxpefac" function will be used.

% Input:   s(ns)      Speech signal
%          fs         Sample frequency (Hz)
%          tinc       Time increment between frames (s) [0.01]
%                     or [start increment end]
%          m          mode
%                     'g' plot graph showing waveform and pitch
%                     'G' plot spectrogram with superimposed pitch using
%                         options pp.sopt [default: 'ilcwpf']
%                     'x' use external files for algorithm parameter
%                         initialization: fxpefac_g and fxpefac_w
%          pp         structure containing algorithm parameters
%
% Outputs: fx(nframe)     Estimated pitch (Hz)
%          tx(nframe)     Time at the centre of each frame (seconds).
%          pv(nframe)     Probability of the frame being voiced
%          fv             structure containing feature vectors
%                           fv.vuvfea(nframe,2) = voiced/unvoiced GMM features

% The output "pv(nframe)" will be used. This is the probability of the frame being voiced.

% Define the filename of the .wav file
filename = 'sp01_white_noise_15dB';
%filename = 'sp01_airport_sn0';

% Read the .wav file
[x, Srate, nbits] = wavread(filename);

% Define the frequency resolution of the initial spectrogram (Hz)
p.fstep=5;              

% Define the maximum frequency of the initial spectrogram (Hz)
p.fmax=4000;            

% Define the bandwidth of the initial spectrogram (Hz)
p.fres = 20;            

% Define the frame increment (s)
p.tinc = 0.01;          

% Define the spectrogram of the mixture
fmin = 0; fstep = p.fstep; 
fmax = p.fmax;

% Frequency resolution (Hz)
fres = p.fres;  

% Use the function "spgrambw" 
% Draw spectrogram, and use: [T,F,B] = spgrambw(s,fs,mode,bw,fmax,db,tinc,ann)
[tx, ~, ~] = spgrambw(x, Srate, fres, [fmin fstep fmax], [], p.tinc);

% Define the number of frames
nframes = length(tx);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Determine if the frame is voiced or unvoiced 

% Use the function "fxpefac": [fx,tx,pv,fv] = fxpefac(s,fs,tinc,m,pp)
[~, ~, pv, ~] = fxpefac(x, Srate);
% We will use "pv(n)"
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Initialize the variables for the for loop
total_noise_ps = [];

% Define the threshold value for the voiced speech probability
threshold = 0.1 * 10^(-4);

% Define the length of a single frame
len = floor(Srate / nframes); 

% Define the window to be used 
win = hamming(len); 

% Define the FFT length
nFFT = 2 * len;

% for loop, for processing all the frames
%for n = 1 : Nframes 
for n = 1 : nframes 
    % We have determined if the frame is voiced or unvoiced
    % We have used the function "fxpefac": [fx,tx,pv,fv] = fxpefac(s,fs,tinc,m,pp)
    % We now use "pv(n)", which is the probability of the frame being voiced
   
    % Perform windowing  
    input_signal = win .* x(n:n+len-1); 

    %Find the Fourier transform (FFT) of a frame
    FFT_input_signal = fft(input_signal, nFFT);    

    % Calculate the magnitude of the FFT
    magnitude_input_signal = abs(FFT_input_signal); 

    % Define the power spectrum of the noise
    power_noise = magnitude_input_signal .^ 2;
    
    % If the frame has only noise
    if (pv(n) < threshold)
        total_noise_ps = [total_noise_ps (power_noise)'];
    
    % If the frame has speech and noise. The frane does not have only noise.
    else
        
        

        
        
        
    end
    
    


    
    


    
%    % Estimate and update the noise power spectrum 
%    if n == 1
%          parameters = initialise_parameters(ns_ps, Srate, method);   
%     else
%         parameters = noise_estimation(ns_ps, method, parameters);
%    end
% 
%    % Find the noise power spectrum
%    noise_ps = parameters.noise_ps;
% 
%    % Update the total noise power spectrum
%    total_noise_ps = [total_noise_ps; noise_ps'];
end










 
 
 

% This method has 7 more frames than the previous method all the times.

% Read the .wav file
[x, Srate, nbits] = wavread(filename);

% Define the frame size (in samples)
% We use a frame length of 320 samples
len = floor(20 * Srate / 1000); 

% We ensure that the frame size is valid
if (rem(len,2) == 1) 
    len = len + 1; 
end

% Define the window overlap (percentage)
PERC = 50; 
len1 = floor(len * PERC / 100);
len2 = len - len1;

% Define the window to be used 
win = hamming(len); 

% Initialize the variables
k = 1;
nFFT = 2 * len;
img = sqrt(-1);
x_old = zeros(len1, 1);
Nframes = floor(length(x) / len2) - 1;
xfinal = zeros(Nframes * len2, 1);










% % Use the function "activlevg" to measure the active speech level 
% % We use: [LEV, AF, FSO] = activlevg(sp,FS,MODE)
% 
% % Inputs: sp     is the speech signal
% %         FS     is the sample frequency in Hz (see also FSO below)
% %         MODE   is a combination of the following:
% %                r - raw omit input filters (default is 200 Hz to 5.5 kHz)
% %                0 - no high pass filter (i.e. include DC)
% %                4 - high pass filter at 40 Hz (but allows mains hum to pass)
% %                1 - use cheybyshev 1 filter
% %                2 - use chebyshev 2 filter (default)
% %                e - use elliptic filter
% %                h - omit low pass filter at 5.5 kHz
% %                d - give outputs in dB rather than power
% %                n - output a normalized speech signal as the first argument
% %                N - output a normalized filtered speech signal as the first argument
% %                l - give additional power level estimates (see below for details)
% %                a - include A-weighting filter
% %                i - include ITU-R-BS.468/ITU-T-J.16 weighting filter
% % 
% % Outputs:
% %     If the "n" option is specified, a speech signal normalized to 0dB will be given as
% %     the first output followed by the other outputs.
% %         LEV    gives the speech level in units of power (or dB if mode='d')
% %                if mode='l' is specified, LEV is a row vector containing:
% %                        [ active-level mean-power mean-noise-power P.56-level harmonic-power-level]
% 
% % Call the function to find the active speech level
% LEV = activlevg(x, Srate, 'dl');
% 
% sa










% Initialize the variables for the for loop
total_noise_ps = [];

% for loop, for processing all the frames
for n = 1 : Nframes 
    % Perform windowing  
    insign = win.*x(k:k+len-1); 

    %Find the Fourier transform (FFT) of a frame
    spec = fft(insign,nFFT);    

    % Calculate the magnitude of the FFT
    sig = abs(spec); 

    % Calculate the power spectrum
    ns_ps = sig.^2;

    % Use the function "activlevg" to measure the active speech level 
    % We use: [LEV, AF, FSO] = activlevg(sp,FS,MODE)

    
    % Define the normalized speech power
    normalized_ns_ps = ns_ps / LEV;
    
    asdfa
    
    % We use the function "histndim" to generate and plot an n-dimensional histogram

    %   Inputs:  X(m,d)   is the input data: each row is one d-dimensiona data point
    %            B(3,d)   specifies the histogram bins.
    %                          B(1,:) gives the number of bins in each dmension [default 10]
    %                          B(2,:) gives the minimum of the first bin in each dimension [default min(X)]
    %                          B(3,:) gives the maximum of the last bin in each dimension [default max(X)]
    %                     If B has only one column, the same values are use for al dimensions 
    %                     If B(1,i)=0 then that dimension will be ignored (and excluded from V)
    %            MODE     is a character string containing a combination of the following:
    %                         'z' for zero base in the 2D plot [default base = min(V)]
    %                         'p' to scale V as probabilities [default actual counts]
    %                         'h' to plot a histogam even if output arguments are present
    % 
    %  Outputs:  V        d-dimensional array containing the histogram counts
    %            T        d-element cell array. d{i} contains the bin boundary values for
    %                     the i'th dimension. The length of d{i} is one more than the number of bins
    %                     in that dimension.
    % 
    %                     Note that if any of B(1,:) are zero then the number of dimensions in V and elements
    %                     of T will be correspondingly reduced.
    % 
    %  Example: histndim(randn(100000,2),[20 -3 3]','pz');
    
    %histndim(normalized_ns_ps, [599902 300 3000], 'p');
    
      
   
   
   
%    % Use the function "gaussmix" to fit a Gaussian mixture pdf to a set of data observations 
%    %     M(k,p)   Mixture means, one row per mixture. 
%    %     V(k,p)   Mixture variances, one row per mixture. 
%    %     W(k,1)   Mixture weights, one per mixture. The weights will sum to unity. 
%       
%    % Create a GMM with k mixtures
%    k = 2; 
%    [m, v, w] = gaussmix(,[],[],k);    
   
   
    
   

   
   
%    % Estimate and update the noise power spectrum 
%    if n == 1
%          parameters = initialise_parameters(ns_ps, Srate, method);   
%     else
%         parameters = noise_estimation(ns_ps, method, parameters);
%    end
% 
%    % Find the noise power spectrum
%    noise_ps = parameters.noise_ps;
% 
%    % Update the total noise power spectrum
%    total_noise_ps = [total_noise_ps; noise_ps'];
end












%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%

% In this part of the code, the "mmse_noise_estimation" function will be used.
% We use: [noise] = mmse_noise_estimation(filename, SPU)

% Define not to use the Speech Presence Probability
SPU = 0;

% Call the function
filename = 'sp02_train_sn5.wav';
noise = mmse_noise_estimation(filename, SPU);

% Define nFFT for the noise signal
nFFT = floor(length(noise)) - 1;

% Find the Fourier transform (FFT) of the noise
spec = fft(noise, nFFT);    

% Calculate the magnitude of the FFT
sig = abs(spec); 

% Calculate the power spectrum
ns_ps = sig.^2;
















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, 2 spectrograms of 2 speech signals will be plotted.
% The function "spgrambw" from the M.Brookes' Voicebox will be used.

% the following function returns the sampled data, the sample rate (Fs) in Hertz and 
% the number of bits per sample (NBITS) used to encode the data in the file.
[input_speech, Fs, NBITS] = wavread('senta_ms.wav');
% The speech signal is: "Betty can see the wagon".

figure; set(gcf,'Color','w');
% Call the function 
spgrambw(input_speech, Fs, 'pJcw');
% Set the title of the image
set(gca,'FontSize',15);
title('Spectrogram of a Speech Signal from a male speaker. The sentence is: "Betty can see the wagon".');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% the following function returns the sampled data, the sample rate (Fs) in Hertz and 
% the number of bits per sample (NBITS) used to encode the data in the file.
[input_speech, Fs, NBITS] = wavread('sentb_ms.wav');
% The speech signal is: "We can't go there today".

figure; set(gcf,'Color','w');
% Call the function 
spgrambw(input_speech, Fs, 'pJcw');
% Set the title of the image
set(gca,'FontSize',15);
title('Spectrogram of a Speech Signal from a female speaker. The sentence is: "We cant go there today".');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, certain noise estimation algorithms will be tested.

% The noise estimation algorithms are the following.
% The algorithms are sorted in chronological order.
% 1) 'mcra2': MCRA-2, Loizou (2006), MCRA = Minima Controlled Recursive averaging,
%    A noise estimation algorithm for highly non-trainary environments. 
% 2) 'conn_freq': Sorensen (2005), Connected time-frequency regions, Speech enhancement with 
%    natural sounding residual noise based on connected time-frequency speech presence regions.

% 3) 'imcra': IMCRA, Cohen (2003), IMCRA = Improved Minima Controlled Recursive averaging,
%    Noise spectrum estimation in adverse enviroments. 
% 4) 'mcra': MCRA algorithm, Cohen (2002), MCRA = Minima Controlled Recursive averaging, 
%    Noise estimation by minima controlled recursive averaging for robust speech enhancement. 

% 5) 'martin': Martin (2001), Minimum tracking algorithm, Noise power spectral density
%    estimation based on optimal smoothing and minimum statistics.
% 6) 'doblinger': Doblinger (1995), Continuous minimal tracking, Computationally efficient 
%    speech enhancement by spectral minima tracking in subbands. 
% 7) 'hirsch': Hirsch (1995), Weighted spectral average, Noise estimation techniques 
%    for robust speech recognition. 

% (1) MCRA-2, Loizou (2006), (use: 'mcra2') 
% Define the method to be used 
method = 'mcra2';
% Define the name of the .wav file
filename = 'sp02_train_sn5.wav';
% Call the function so as to find the noise power spectrum for each frame
[mcra2_noise_ps] = noise_parameters(filename, method);

% (2) Connected time-frequency regions, Sorensen (2005), (use: 'conn_freq')  
% Define the method to be used 
method = 'conn_freq';
% Define the name of the .wav file
filename = 'sp02_train_sn5.wav';
% Call the function so as to find the noise power spectrum for each frame
[conn_freq_noise_ps] = noise_parameters(filename, method);

% (3) IMCRA, Cohen (2003), (use: 'imcra')
% Define the method to be used 
method = 'imcra';
% Define the name of the .wav file
filename = 'sp02_train_sn5.wav';
% Call the function so as to find the noise power spectrum for each frame
[imcra_noise_ps] = noise_parameters(filename, method);

% (4) MCRA, Cohen (2002), (use: 'mcra') 
% Define the method to be used 
method = 'mcra';
% Define the name of the .wav file
filename = 'sp02_train_sn5.wav';
% Call the function so as to find the noise power spectrum for each frame
[mcra_noise_ps] = noise_parameters(filename, method);

% (5) Minimum tracking algorithm, Martin (2001), (use: 'martin')
% Define the method to be used 
method = 'martin';
% Define the name of the .wav file
filename = 'sp02_train_sn5.wav';
% Call the function so as to find the noise power spectrum for each frame
[martin_noise_ps] = noise_parameters(filename, method);

% (6) Continuous minimal tracking, Doblinger (1995), (use: 'doblinger')  
% Define the method to be used 
method = 'doblinger';
% Define the name of the .wav file
filename = 'sp02_train_sn5.wav';
% Call the function so as to find the noise power spectrum for each frame
[doblinger_noise_ps] = noise_parameters(filename, method);

% (7) Weighted spectral average, Hirsch (1995), (use: 'hirsch') 
% Define the method to be used 
method = 'hirsch';
% Define the name of the .wav file
filename = 'sp02_train_sn5.wav';
% Call the function so as to find the noise power spectrum for each frame
[hirsch_noise_ps] = noise_parameters(filename, method);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, the noise will be estimated with a Voice Activity Detector (VAD)

% Define the file with the noisy speech signal
filename = 'sp02_train_sn5.wav';

% The function "VAD_noise_estimation" is called
noise_estimation = VAD_noise_estimation(filename);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% Paper: "Unbiased MMSE-Based Noise Power Estimation with Low Complexity and Low Tracking Delay"
% Author: T.Gerkmann and R.C.Hendriks

% Define the filename of the .wav file
filename = 'sp01_white_noise_15dB';
%filename = 'sp01_airport_sn0';

% Read the .wav file
[x, Srate, ~] = wavread(filename);

% Call the "noisePowProposed" function of T.Gerkmann and R.C.Hendriks
noisePowMat = noisePowProposed(x, Srate);
% noisePowMat has dimensions: (noise psd) x (frames)
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% Paper: "MMSE BASED NOISE PSD TRACKING WITH LOW COMPLEXITY"
% Author: T.Gerkmann and R.C.Hendriks

% Add datapath
addpath ./TabGenGam

% Define the filename of the .wav file
filename = 'sp01_white_noise_15dB';
%filename = 'sp01_airport_sn0';

% Read the .wav file
[x, Srate, ~] = wavread(filename);

% Call the "noisePowProposed" function of T.Gerkmann and R.C.Hendriks
[shat, noise_psd_matrix, T] = noise_psd_tracker(x, Srate);
% noisePowMat has dimensions: (noise psd) x (frames)
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, the algorithms are integrated into a spectral subtraction algorithm.

% The function "specsub_ns" is used: specsub_ns(filename,method,outfile)
%         method - noise estimation algorithm:
%                  'martin'    = Martin''s minimum tracking algorithm
%                  'mcra'      = Minimum controlled recursive average algorithm (Cohen) 
%                  'mcra2'     = Variant of Minimum controlled recursive average algorithm (Rangachari, Loizou)
%                  'imcra'     = Improved Minimum controlled recursive average algorithm (Cohen)
%                  'doblinger' = Continuous spectral minimum tracking(Doblinger) 
%                  'hirsch'    = Weighted spectral average (Hirsch, Ehrilcher) 
%                  'conn_freq' = Connected frequency regions (Sorensen, Andersen)

% Define the input parameters
filename = 'sp01_white_noise_15dB.wav';
outfile = 'sp01_white_noise_15dB_spectralSub_martin.wav';

% Call the spectral subtraction algorithm
% We use the Martin''s minimum tracking algorithm
specsub_ns(filename, 'martin', outfile);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use PESQ: scores = pesq(cleanFile, enhancedFile)

%         cleanFile     - clean input file in .wav format sampled at 
%                         sampling frequency Fs=8 kHz or Fs=16 kHz
%                         for narrowband or wideband assessment,
%                         respectively.
%
%         enhancedFile  - enhanced output file in .wav format sampled
%                         at same sampling frequency as the cleanFile
%
%         scores        - For narrowband speech, two scores are returned,
%                         one for the raw PESQ value [1] (first value) and 
%                         one for the MOS-mapped score value [2] (second value).
%                         For wideband speech, only the MOS-mapped value
%                         is returned [3].

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';

% Change the sampling frequency of the clean speech file
[x, Fs, nbits] = wavread('sp01.wav');
[P, Q] = rat(16000/Fs);
x_new = resample(x,P,Q);
wavwrite(x_new, 16000, nbits, 'sp01_16kHz.wav');

% Define the clean speech file
cleanFile = 'sp01_16kHz.wav';

% We call the "PESQ" function
scores = pesq(cleanFile, enhancedFile);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the Cepstrum Distance Objective Speech Quality Measure.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the "comp_cep" function
cep_mean = comp_cep(cleanFile, enhancedFile);
% The "cep_mean" is the computed cepstrum distance measure.
% The cepstrum measure is limited in the range [0, 10].
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the frequency weighted SNRseg Objective Speech Quality Measure.

% We use the function "comp_fwseg".
% The function implements the frequency-weighted SNRseg measure 
% using the clean spectrum as a weighting function.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the "comp_fwseg" function
fwseg_dist = comp_fwseg(cleanFile, enhancedFile);
% The "fwseg_dist" is the computed frequency weighted SNRseg in dB
% The larger the "fwseg_dist" is, the better the spectral subtraction algorithm is.
% Hence, the better the noise estimation algorithm is.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the MARS Frequency-variant fwSNRseg objective speech quality measure.

% We use the function "comp_fwseg_mars".
% The function implements the frequency-variant fwSNRseg measure based on MARS analysis.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the function to perform the MARS frequency-variant fwSNRseg
[sig, bak, ovl] = comp_fwseg_mars(cleanFile, enhancedFile);
%         sig           - predicted rating [1-5] of speech distortion
%         bak           - predicted rating [1-5] of noise distortion
%         ovl           - predicted rating [1-5] of overall quality
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the frequency-variant fwSNRseg Objective Speech Quality Measure.

% The function "comp_fwseg_variant" will be used.
% The function implements the frequency-variant fwSNRseg measure.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the function to perform the frequency-variant fwSNRseg.
[sig, bak, ovl] = comp_fwseg_variant(cleanFile, enhancedFile);
%         sig           - predicted rating [1-5] of speech distortion
%         bak           - predicted rating [1-5] of noise distortion
%         ovl           - predicted rating [1-5] of overall quality
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the Itakura-Saito (IS) Objective Speech Quality Measure

% The function "comp_is" will be used.
% The function implements the Itakura-Saito distance measure.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the function to find the IS.
is_mean = comp_is(cleanFile, enhancedFile);
%         IS            - computed Itakura Saito measure
%         The IS measure is within the range [0, 100].
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the Log Likelihood Ratio (LLR) Objective Speech Quality Measure.

% The function "comp_llr" will be used.
% The function implements the Log Likelihood Ratio Measure.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the function to find the LLR.
llr_mean = comp_llr(cleanFile, enhancedFile);
%         llr           - computed likelihood ratio
%         The LLR measure is within the range [0, 2].
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the Segmental Signal-to-Noise Ratio Objective Speech Quality Measure.

% The function "comp_snr" will be used.
% The function implements the segmental signal-to-noise ratio.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the function to find the Segmental SNR (SegSNR).
[snr_mean, segsnr_mean] = comp_snr(cleanFile, enhancedFile);
%         SNRovl        - overall SNR (dB)
%         SNRseg        - segmental SNR (dB)

% The second value is the segmental signal-to-noise ratio (1 seg-snr per frame of input).  
% The segmental SNR is within the range 35 dB and -10 dB.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the Weighted Spectral Slope (WSS) Objective Speech Quality Measure.

% The function "comp_wss" will be used.
% The function implements the WSS.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the function to find the Segmental SNR (SegSNR).
wss_dist = comp_wss(cleanFile, enhancedFile);
% The "wss_dist" is the computed spectral slope distance.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We evaluate the spectral subtraction algorithm.
% We use the Composite Objective Speech Quality Measure

% The function "composite" will be used.
% The function implements the Composite measure.

% Define the input parameters
enhancedFile = 'sp01_white_noise_15dB_spectralSub_martin.wav';
cleanFile = 'sp01_16kHz.wav';

% We call the function to find the Composite measure.
[Csig, Cbak, Covl] = composite(cleanFile, enhancedFile);
%         Csig           - predicted rating [1-5] of speech distortion
%         Cbak           - predicted rating [1-5] of noise distortion
%         Covl           - predicted rating [1-5] of overall quality
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% Create modulated Gaussian noise

% This is based on T.Gerkmann and R.C.Hendriks (2012) paper.
modulated_white_Gaussian_noise = 1 + 0.5*sin(2*pi*(0.5/16000)*(1:16000));
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, we will test the "MatlabADT" functionality.

% Add the path to the "MatlabADT" programme
addpath ./MatlabADT
db = ADT('timit','E:\project_draft','setup');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% This call to the query function will return the wave data of the first 30 words 'she'
% or 'it' form dialect 'dr1' in the form of a cell array
[wave fs metadata] = query(db,'dialect','dr1','word',{'she','it'},30);

% Define one word
oneWord = wave{1};
%soundsc(oneWord,16000)

% plot the result for one word
figure; set(gcf, 'Color', 'w');
plot(oneWord);
set(gca, 'FontSize', 15);
title('Time Domain Representation of the word SHE.');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% This call to the query function will return the SA1 sentences
[SA1_wave2 fs2 metadata2] = query(db, 'sentence', 'SA1', 30);
% we use: SA1_wave2{1}

% This call to the query function will return the SA2 sentences
[SA2_wave2 fs2 metadata2] = query(db, 'sentence', 'SA2', 30);
% we use: SA2_wave2{1}
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------












%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% Create a GUI table for R.Martin's MS algorithm

% We use R.Martin's MS algorithm
method = 'martin';

filename = 'sp04_babble_sn10';
outfile = strcat('output_', filename);

noise_ps_martin = specsub_ns(filename, method, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp04';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = pesq_val(2);

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp04_babble_sn5';

outfile = strcat('output_', filename);

%method = 'mcra';
method = 'martin';

noise_ps_martin = specsub_ns(filename, method, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp04';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = pesq_val(2);

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final2 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp04_babble_sn15';

outfile = strcat('output_', filename);

method = 'martin';

noise_ps_martin = specsub_ns(filename, method, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp04';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = pesq_val(2);

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final3 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp04_babble_sn0';

outfile = strcat('output_', filename);

method = 'martin';

noise_ps_martin = specsub_ns(filename, method, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp04';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = pesq_val(2);

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final4 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

%create a table with the values calculated
T = table(final4, final2, final3, final);
T.Properties.VariableNames = {'TestResult', 'TestResult2', 'TestResult3', 'TestResult4'};

%print the results
disp(['The following table shows the results that we just calculated.']);
disp(T);
disp(['   ']); disp(['   ']);

%set the figure data
f = figure; set(gcf,'Color','w');
set(f, 'Position', [500 400 300 100]);
dat =  {T{1,1}, T{1,2}, T{1,3}, T{1,4};...
        T{2,1}, T{2,2}, T{2,3}, T{2,4};...   
        T{3,1}, T{3,2}, T{3,3}, T{3,4};...
        T{4,1}, T{4,2}, T{4,3}, T{4,4};...
        T{5,1}, T{5,2}, T{5,3}, T{5,4};...
        T{6,1}, T{6,2}, T{6,3}, T{6,4};...
        T{7,1}, T{7,2}, T{7,3}, T{7,4};...
        T{8,1}, T{8,2}, T{8,3}, T{8,4};...
        T{9,1}, T{9,2}, T{9,3}, T{9,4};...
        T{10,1}, T{10,2}, T{10,3}, T{10,4};...
        T{11,1}, T{11,2}, T{11,3}, T{11,4};...
        T{12,1}, T{12,2}, T{12,3}, T{12,4};...
        T{13,1}, T{13,2}, T{13,3}, T{13,4};...
        T{14,1}, T{14,2}, T{14,3}, T{14,4};...
        T{15,1}, T{15,2}, T{15,3}, T{15,4};...
        T{16,1}, T{16,2}, T{16,3}, T{16,4};...
        T{17,1}, T{17,2}, T{17,3}, T{17,4};...
        T{18,1}, T{18,2}, T{18,3}, T{18,4};};

%set the table as a figure
columnname =   {'Results for 15 dB SNR', 'Results for 10 dB SNR', 'Results for 5 dB SNR', 'Results for 0 dB SNR'};
rnames = {'SNR mean', 'Segmental SNR', 'WSS mean', 'LLR mean', 'IS mean', 'CEP mean', ...
    'fw SNR segmental', 'SIG' 'BAK', 'OVL', 'SIG2' 'BAK2', 'OVL2', 'raw PESQ', 'MOS-valued PESQ', 'Csig', 'Cbak', 'Covl'};
columnformat = {'numeric'}; 
t = uitable('Units','normalized','Position',...
            [0.05 0.05 0.755 0.87], 'Data', dat,... 
            'ColumnName', columnname,...
            'RowName', rnames, ...
            'ColumnFormat', columnformat); 
        
%set the title
%set(gca,'FontSize',15); 
title('Table: R.Martin MS Algorithm. Only babble noise is used. Only one speaker is tested.');        
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------











%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% In this part of the code, the noisy signals are created.

% % The "NOIZEUS_Clean_Speech_Database" for the clean signal (sp01.wav) is used.
% % The "Noise_Database" for the noise signal (white_noise.wav) is used.
% % The noisy output signal has a SNR of 0 dB.
% addnoise_asl('sp01.wav', 'white_noise.wav', 'sp01_white_noise_15dB.wav', 15);
% 
% % The "NOIZEUS_Clean_Speech_Database" for the clean signal (sp01.wav) is used.
% % The "Noise_Database" for the noise signal (white_noise.wav) is used.
% % The noisy output signal has a SNR of 5 dB.
% addnoise_asl('sp01.wav', 'white_noise.wav', 'sp01_white_noise_10dB.wav', 10);
% 
% % The "NOIZEUS_Clean_Speech_Database" for the clean signal (sp01.wav) is used.
% % The "Noise_Database" for the noise signal (white_noise.wav) is used.
% % The noisy output signal has a SNR of 5 dB.
% addnoise_asl('sp01.wav', 'white_noise.wav', 'sp01_white_noise_5dB.wav', 5);
% 
% % The "NOIZEUS_Clean_Speech_Database" for the clean signal (sp01.wav) is used.
% % The "Noise_Database" for the noise signal (white_noise.wav) is used.
% % The noisy output signal has a SNR of 0 dB.
% addnoise_asl('sp01.wav', 'white_noise.wav', 'sp01_white_noise_0dB.wav', 0);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % We need to use the filename and change it every time
% filename = 'ns_sps_0dB_S_01_01';
% 
% % Outter for loop
% for i = 1 : 10
%     % Inner for loop
%     for k = 1 : 72
%         % We change the filename based on i
%         % If statement for changing the filename based on i
%         if (i ~= 10)
%             filename(18) = num2str(i);
%         else
%             filename(17) = num2str(1);
%             filename(18) = num2str(0);
%         end
% 
%         % We change the filename based on k
%         % We use a temporary variable
%         temp_variable = num2str(k);
%     
%         % We now change the filename based on k
%         if (k < 10) 
%             filename(14) = num2str(0);
%             filename(15) = num2str(temp_variable(1));
%         else
%             filename(14) = num2str(temp_variable(1));
%             filename(15) = num2str(temp_variable(2));
%         end
%         
%         % Read the .wav file
%         [x, Fs, nbits] = wavread(filename);
%         
%         % Resample the signal
%         [P, Q] = rat(16000/Fs);
%         x_new = resample(x,P,Q);
%         
%         % Write the result
%         wavwrite(x_new, 16000, nbits, filename);
%     end
% end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% We assume that we know the distribution of the noisy speech signal
% We assume that it is Gaussian 

% We define the parameters
x = [-110 : 0.01 : -40];
mn = -90;
var = 2;

% We assume that the pdf is Gaussian 
norm = normpdf(x, mn, var);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Plot the pdf
figure; set(gcf,'Color','w');
plot(x,norm);

% Set the title of the image
set(gca,'FontSize',15);
title('Normal pdf for the Noisy Speech Signal.');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

norm = exp(norm);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------


% % Create a GMM model
% mu = [1 2;-3 -5];
% sigma = cat(3,[2 0;0 .5],[1 0;0 1]);
% p = ones(1,2)/2;
% obj = gmdistribution(mu,sigma,p);
% 
% %ezsurf(@(x,y)pdf(obj,[x y]),[-10 10],[-10 10])
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% x = -10:0.1:10;
% y = -10:0.1:10;
% % plot(pdf(obj,[x' y']))
% 
% norm = pdf(obj,[x' y']);
% norm = exp(norm);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Find the autocorrelation R matrix

% In order to find a bound on |t|, we plot the estimate for various |t|
% We recommend that |t| is at about 20-30% of N.
% We set the upper limit to 30% of N.

% Define N
N = length(norm);

% Define x
x = norm;
%x = norm';

% Compute the estimate of the ACF with |t| = (0.25*N)
acf_temp = 0; acfSizeMatrix = 0; acf_tempTotal = 0; acf_total = 0;
t = 0.25 * N; 
acf_temp = acfZoom2(x, t);

% We use the fact that the ACF is symmetric on the y-axis.
acfSizeMatrix = size(acf_temp);

% Find the samples from lag 1 to lag acfSizeMatrix(1,2)
acf_tempTotal = acf_temp(1, 2:acfSizeMatrix(1,2));

% Invert the samples
for i = 1 : (acfSizeMatrix(1,2)-1)
    acf_total(1, i) = acf_tempTotal(1, ((acfSizeMatrix(1,2)-1)-i+1));
end

% Find the acf estimate when zoom in occurs
acf_total = [acf_total acf_temp];
acfSizeMatrix = size(acf_total);

% Plot the results
figure; set(gcf,'Color','w');
stem(((-(acfSizeMatrix(1,2)+1)/2)+1) : (((acfSizeMatrix(1,2)+1)/2)-1), acf_total);

% Name the x and y axes
set(gca,'FontSize',15); 
xlabel('Lag'); ylabel('ACF estimate');
title('ACF estimate. We use the zoom factor |t|. We suggest |t| upper limit to be 20-30% of N.');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use R(0) = acf_total(((acfSizeMatrix(1,2)+1)/2)-1)
% We use R(0) = acf_temp(1)
R_matrix = acf_temp;

% define n
n = 1;

% find the coefficient c_1
c_1 = inv(R_matrix(1)) * R_matrix(2);

% mean square error
temp = conj(R_matrix(2)') * inv(R_matrix(1)) * R_matrix(2);
sigma_squared = R_matrix(1) - temp;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% n = 2

% initilize the matrices
r_WienerHopf = [];
R_matrix_WienerHopf = [];

% define n
n = 2;

for i = 1 : n
    % define the Bessel function
    r_WienerHopf = [r_WienerHopf; R_matrix(1+n)];
end

for i = 1 : n
    for j = 1 : n
        % define the R matrix
        R_matrix_WienerHopf(i,j) = R_matrix(1+abs(i-j));
    end
end

% find the coefficients c
c = inv(R_matrix_WienerHopf) * r_WienerHopf;

% mean square error
temp = (r_WienerHopf') * inv(R_matrix_WienerHopf) * r_WienerHopf;
sigma_squared =  R_matrix(1) - temp;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% for all n

% initilize the matrices
total_sigma_squared = [];

% for loop, define n
for n = 1 : 20
    % initilize the matrices
    r_WienerHopf = [];
    R_matrix_WienerHopf = [];

    for i = 1 : n
        % define the r_WienerHopf function
        r_WienerHopf = [r_WienerHopf; R_matrix(1+n)];
    end

    for i = 1 : n
        for j = 1 : n
            % define the R matrix
            R_matrix_WienerHopf(i,j) = R_matrix(1+abs(i-j));
        end
    end

    % find the coefficients c
    c = inv(R_matrix_WienerHopf) * r_WienerHopf;

    % mean square error
    temp = (r_WienerHopf') * inv(R_matrix_WienerHopf) * r_WienerHopf;
    sigma_squared =  R_matrix(1) - temp;
    
    % store the result
    total_sigma_squared = [total_sigma_squared sigma_squared];
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% plot the results
figure; set(gcf,'Color','w');
stem(total_sigma_squared);
xlim([0.5 20.5]);

% set the title
set(gca,'FontSize',15);
title('Minimum Mean Square Error');
ylabel('MMSE, sigma_n^2'); xlabel('n');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% plot the results
figure; set(gcf,'Color','w');
plot(total_sigma_squared);
xlim([0.5 20.5]);

% set the title
set(gca,'FontSize',15);
title('Minimum Mean Square Error');
ylabel('MMSE, sigma_n^2'); xlabel('n');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------













%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

filename = 'sp01_airport_sn0';
[norm, Fs] = readwav(filename);

[filtdata, f_filtdata, ~] = wienerFilter(norm,norm,1,0,Fs); 


observation = norm;
ideal = filtdata;

% estimate noise from ideal
noise = observation-ideal;

% work out how long to make FFT
N=length(observation);

% Wiener filter
Sf2=real(fft(ideal,N*2-1)).^2;   % Smeared ideal
Nf2=real(fft(noise,N*2-1)).^2;   % noise
Cf=real(fft(observation,N*2-1)); % ~= sqrt(Sf2+Nf2); % Corrupted ideal

H=Sf2./(Sf2+Nf2);              % Optimal filter

Yhat=H.*Cf;                  % 'uncorrupted' ideal estimate ...

yhat=real(ifft(Yhat));           % ..... in time domain

% ...compensate for FFT being two sided in matlab   
yhat=yhat(1:length(observation)); 


for i = 1 : 3
    observation = norm;
    ideal = yhat;

    % estimate noise from ideal
    noise = observation-ideal;

    % work out how long to make FFT
    N=length(observation);

    % Wiener filter
    Sf2=real(fft(ideal,N*2-1)).^2;   % Smeared ideal
    Nf2=real(fft(noise,N*2-1)).^2;   % noise
    Cf=real(fft(observation,N*2-1)); % ~= sqrt(Sf2+Nf2); % Corrupted ideal

    H=Sf2./(Sf2+Nf2);              % Optimal filter

    Yhat=H.*Cf;                  % 'uncorrupted' ideal estimate ...

    yhat=real(ifft(Yhat));           % ..... in time domain

    % ...compensate for FFT being two sided in matlab   
    yhat=yhat(1:length(observation)); 
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------








% (3) IMCRA, Cohen (2003), (use: 'imcra')
% Define the method to be used 
method = 'imcra';
% Define the name of the .wav file
filename = 'sp01_train_sn10';
% Call the function so as to find the noise power spectrum for each frame
[imcra_noise_ps] = noise_parameters(filename, method);






% Define the file name
filename = 'sp01_train_sn10';

% We initialize the matrices
total_power_1 = [];


a = 0.1;
%a = 0.2; 
%a = 0.3; 
%a = 0.4; 
%a = 0.5;
%a = 0.6;
%a = 0.7;
%a = 0.8 
%a = 0.9;

% Outter for loop
for outter_loop = 1 : 1
    for inner_loop = 1 : 1
%         % We change the filename based on inner_loop
%         % We use a temporary variable
%         temp_variable = num2str(inner_loop);
% 
%         % We now change the filename based on inner_loop
%         if (inner_loop < 10) 
%             filename(3) = num2str(0);
%             filename(4) = num2str(temp_variable(1));
%         else
%             filename(3) = num2str(temp_variable(1));
%             filename(4) = num2str(temp_variable(2));
%         end
% 
%         % If statement to change the filename
%         if (outter_loop == 1) 
%             filename(14) = num2str(0);
%         elseif (outter_loop == 2) 
%             filename(14) = num2str(5);
%         elseif (outter_loop == 3) 
%             filename(14) = num2str(1);
%             filename(15) = num2str(0);
%         else 
%             filename(14) = num2str(1);
%             filename(15) = num2str(5);
%         end

        % Read the .wav file
        [x, Srate, bits] = wavread(filename);	

        % Define the frame size in samples
        len=floor(20*Srate/1000); 
        if rem(len,2)==1, len=len+1; end;

        % window overlap in percent of frame size
        PERC=50; 
        len1=floor(len*PERC/100);
        len2=len-len1;

        % define window
        win=hanning(len);  

        % normalize window for equal level output 
        win = win*len2/sum(win);  

        % Noise magnitude calculations - assuming that the first 6 frames is noise/silence
        nFFT=2*len;
        j=1;
        noise_mean=zeros(nFFT,1);
        for k=1:6
            noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
            j=j+len;
        end
        noise_mu=noise_mean/6;
        noise_mu2=noise_mu.^2;

        %--- allocate memory and initialize various variables
        k=1;
        img=sqrt(-1);
        x_old=zeros(len1,1);
        Nframes=floor(length(x)/len2)-1;
        xfinal=zeros(Nframes*len2,1);

        % --------------- Initialize parameters ------------
        k=1;
        aa=0.98;
        eta= 0.15;
        mu=0.98;
        c=sqrt(pi)/2;
        qk=0.3;
        qkr=(1-qk)/qk;
        ksi_min=10^(-25/10); 

        % Main for loop
        for n = 1 : Nframes

            insign=win.*x(k:k+len-1);

            %--- Take fourier transform of  frame
            spec=fft(insign,nFFT);
            sig=abs(spec); % compute the magnitude
            sig2=sig.^2;
            
            power_frame = sig2;
            
            if (n == 1)
                first_hypothesis = 0.5;
                second_hypothesis = 0.5;

                sigma_1 = 65;
                sigma_2 = 38;
                
                %a = 0.3;
                
                p_1 = 0.5;
                p_2 = 0.5;

                previous_1 = power_frame(1);
                previous_2 = power_frame(1);
            end

            for i = 1 : 640

                y_squared = power_frame(i);
                energy_frame = sqrt(y_squared);
                
                sigma_1 = previous_1;
                sigma_1 = previous_2;

                inv_1 = (1/sigma_1) * (exp(- (y_squared / sigma_1)));
                inv_2 = (1/sigma_2) * (exp(- (y_squared / sigma_2)));

                % We use 2 hypothesis
                first_hypothesis = (first_hypothesis * inv_1) / ((first_hypothesis * inv_1) + (second_hypothesis * inv_2));

                second_hypothesis = (second_hypothesis * inv_2) / ((first_hypothesis * inv_1) + (second_hypothesis * inv_2));

                expected_1 = first_hypothesis * y_squared + (1 - first_hypothesis) * previous_1;

                expected_2 = second_hypothesis * y_squared + (1 - second_hypothesis) * previous_2;

                if (expected_1 < 0)
                    expected_1 = 0;
                elseif (expected_1 > energy_frame)
                    expected_1 = energy_frame;
                end

                if (expected_2 > 1000 * energy_frame)
                    expected_2 = 1000 * energy_frame;
                elseif (expected_2 < 2 * energy_frame)
                    expected_2 = 2 * energy_frame;
                end

                p_2 = a * p_2 * second_hypothesis + (1 - a * p_2) * p_2;

                if (p_2 < 0)
                    p_2 = 0;
                elseif (p_2 > 0.2)
                    p_2 = 0.2;
                end

                p_1 = 1 - p_2;

                previous_1 = a * p_1 * expected_1 + (1 - a * p_1) * previous_1;
                previous_2 = a * p_2 * expected_2 + (1 - a * p_2) * previous_2;
            end

            % Store the variables
            total_power_1(n,i) = previous_1 * p_1 + previous_2 * p_2;
        end
    end
end







x_algorithm = [];

for i = 1 : size(total_power_1,1)
    x_algorithm = [x_algorithm total_power_1(i,:)];
end

x_reference = [];

for i = 1 : size(imcra_noise_ps,1)
    x_reference = [x_reference imcra_noise_ps(i,:)];
end

%RMS error: for loop
RMS_error_1 = 0;
n = length(x_reference);

%for loop
for i = 1 : n
    RMS_error_1 = RMS_error_1 + ((x_reference(i) - x_algorithm(i))^2);
end

RMS_error_1 = sqrt((1/n) * RMS_error_1);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% % We need to use the filename and change it every time
% filename = 'clean_S_01_01';
% 
% % Initialize the variables for the for loop
% store11 = [];
% 
% % Read the .wav file
% [x, Srate, bits] = wavread( filename);	
% 
% % Initialize the matrix
% total_noise_power_frame = [];
% 
% % Initialize the matrix
% total_speech_active_frames_noisy = [];
% 
% % Define fs
% Srate = 16000;
% 
% % Define the frame size in samples
% % len = txinc_samples*10;
% len = 320;
% 
% f=enframe(x,320, 320/2);
% 
% % window overlap in percent of frame size
% %PERC = 90;
% PERC = 50;
% 
% len1 = 320/2;
% len2 = 320/2;
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Define the parameters
% Nframes = floor(length(x)/len2)-1;
% nFFT = 2 * len;
% 
% x_old = zeros(len1,1);
% xfinal = zeros(Nframes*len2, 1);
% 
% noise_only_final = zeros(Nframes*len2, 1);
% n_old = zeros(len1,1);
% 
% k = 1;
% img = sqrt(-1);
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% store_smoothed_p = [];
% transient_presence_prob = 0.3;
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Define window
% win = hamming(length(f(1,:)'));
% 
% %--- allocate memory and initialize various variables
% k=1;
% img=sqrt(-1);
% x_old=zeros(len1,1);
% Nframes=floor(length(x)/len2)-1;
% xfinal=zeros(Nframes*len2,1);
% 
% % --------------- Initialize parameters ------------
% k=1;
% aa=0.98;
% eta= 0.15;
% mu=0.98;
% c=sqrt(pi)/2;
% qk=0.3;
% qkr=(1-qk)/qk;
% ksi_min=10^(-25/10); 
% 
% power_noise = [];
% index_power = [];
% index_power2 = [];
% 
%         % Main for loop
%         for n = 1 : Nframes
%             % Define the input signal
%             insign = f(n,:)' .* win;
% 
%             % Take Fourier transform of  frame
%             spec = fft(insign, nFFT);
%         
%             sig = abs(spec);
%             
%             sig2 = sig.^2;
% 
%             if (n == 1)
%                noise_mu2 = sig2;
%             end
%             
%             % posteriori SNR
%             gammak=min(sig2./noise_mu2,40);  
% 
%             if n==1
%                 ksi=aa+(1-aa)*max(gammak-1,0);
%             else
%                 ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);   
% 
%                 % decision-direct estimate of a priori SNR
%                 ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%             end
% 
%             log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%             vad_decision = sum(log_sigma_k)/nFFT;    
% 
%             % Noise only frame found
%             if (vad_decision < eta) 
%                 % store power 
%                 power_noise = [power_noise; spec'];
% 
%                 % store the index
%                 index_power2 = [index_power2 n];
%                 
%             else
%                 store11 = [store11; spec'];
%                 
%             end
% 
%             vk=ksi.*gammak./(1+ksi);
%             j0 = besseli(0,vk/2);
%             j1 = besseli(1,vk/2);
% 
%             C=exp(-0.5*vk);
%             A=((c*(vk.^0.5)).*C)./gammak;
%             B=(1+vk).*j0+vk.*j1;
%             hw=A.*B;
% 
%             sig=sig.*hw;
% 
%             % save for estimation of a priori SNR in next frame
%             Xk_prev=sig.^2;  
% 
%             xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%             xi_w= real( xi_w);
% 
%             xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%             x_old= xi_w(len1+ 1: len);
% 
%             k=k+len2;
% 
%             
%         end
% 
% 
% store11 = real(store11);
% 
% % we use: store11
% [afa fasfa] = size(store11);
% variance_hat = [];
% weight_hat = [];
% 
% for i = 1 : fasfa
%     % the variance is the same in the real and imaginary parts    
%     data = store11(:,i);
% 
%     data_mainData = data;
%     clear data;
% 
%     nbins = 100;
%     [y,x] = hist(data_mainData, nbins);
%     [sigma, mu] = gaussfit_new(x, y);
% 
% %     data_mainData = data_mainData - mu;
% %     [y,x] = hist(data_mainData, nbins);
% %     [sigma, mu] = gaussfit_new(x, y);
% 
%     % in the end, we use "sigma" only
%     clear mu;
%     clear y;
%     clear x;
%     clear data;
%     clear data_mainData;
% 
%     % we find variance 
%     sigma_squared = sigma^2;
% 
%     % store the variance
%     variance_hat = [variance_hat; sigma_squared];
% end
% 
% % the variance is the same in the real and imaginary parts    
% % we find the mean power parameter
% variance_hat = 2 * variance_hat;
% 
% % the 3 parameters are: variance_hat, weight_hat
% 
% % we find the mean power parameter
% new_mean_power_parameter2 = variance_hat';
% %new_mean_power_parameter2 = [new_mean_power_parameter2 flipud(new_mean_power_parameter2(2:end-1))];
% 
% % Store the variable
% % save('new_clean_speech_stored_data.mat', 'new_mean_power_parameter2');
% save('final_main_new_new_clean_speech_stored_data.mat', 'new_mean_power_parameter2');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%

% filename = 'sp01_airport_sn5';
% 
% %total_noise_ps = noise_parameters_my_proposed_algorithm( filename );
% total_noise_ps = new_noise_parameters_my_proposed_algorithm( filename );
% 
% true_total_noise_ps = true_noise_parameters_my_proposed_algorithm( filename, 'sp01_16kHz' );
% 
% [MSE_3, RMS_error_3] = MSE_function(true_total_noise_ps, total_noise_ps);
% 
% % find non-relative MSE 
% % find the MedSE
% [MSE_notRelative_3, MedSE_final_3] = MedSE_function(true_total_noise_ps, total_noise_ps);
% 
% 
% safsafasfa
% 
% 
% 
% 
% % find MSE
% 
% MSE_3_total = [];
% MSE_notRel_3_total = [];
% MSE_med_3_total = [];
% 
% filename = 'SA1.WAV';
% 
% % Read the noisy speech file
% [input_main_noisy_speech_signal_z, fs] = readsph((filename));
% 
% filename99 = 'buccaneer1';
% 
% % Read the noise file
% [true, fs2] = readwav((filename99));
% 
% snr = 5;
% 
% % Add noise from column vector n() at sample frequency fn with random start
% % Sample and wrapping around as required with a vorbis cross-fade
% z = v_addnoise(input_main_noisy_speech_signal_z,fs,snr,'',true,fs2); 
% 
% z434 = v_addnoise(input_main_noisy_speech_signal_z,fs,snr,'x',true,fs2); 
% 
% input_noise_noiseOnly_signal = z434(:,2);
% 
% % R. Martin MS
% % MS algorithm 
% [MSE_12, MSE_notRelative_12, MedSE_final_12] = main_proposed_algorithm_project( 8, z, input_noise_noiseOnly_signal );
% 
% 
% 
% 
% 
% MSE_3_total = [MSE_3_total MSE_12];
% MSE_notRel_3_total = [MSE_notRel_3_total MSE_notRelative_12];
% MSE_med_3_total = [MSE_med_3_total MedSE_final_12];   
% 
% MSE_3_final_total = mean(MSE_3_total);
% final_MSE_notRel_3_total = mean(MSE_notRel_3_total);
% final_MSE_med_3_total = mean(MSE_med_3_total);
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------

% find MSE

% We now use the "NOIZEUS_Noisy_Speech_Signal_Database_16kHz" database.
%filename = 'sp01_airport_sn0';
%filename = 'sp01_airport_sn5';

filename = 'sp01_exhibition_sn5';
%filename = 'sp01_airport_sn10';
%filename = 'sp01_airport_sn15';

MSE_3_total = [];

% Outter for loop
%for outter_loop = 1 : 4
for outter_loop = 1 : 4
    %for inner_loop = 1 : 30
    for inner_loop = 1 : 30
        
%         if (inner_loop == 25)
%            inner_loop = 26; 
%         end
        
        if (outter_loop == 3 && inner_loop == 4)
            inner_loop = 5;
        end


        % We change the filename based on inner_loop
        % We use a temporary variable
        temp_variable = num2str(inner_loop);

        % We now change the filename based on inner_loop
        if (inner_loop < 10) 
            filename(3) = num2str(0);
            filename(4) = num2str(temp_variable(1));
        else
            filename(3) = num2str(temp_variable(1));
            filename(4) = num2str(temp_variable(2));
        end

        % If statement to change the filename
        if (outter_loop == 1) 
            %filename(16) = num2str(0);
            filename(19) = num2str(0);
        elseif (outter_loop == 2) 
            %filename(16) = num2str(5);
            filename(19) = num2str(5);
        elseif (outter_loop == 3) 
            %filename(16) = num2str(1);
            filename(19) = num2str(1);
            
            %filename(17) = num2str(0);
            filename(20) = num2str(0);
        else 
            %filename(16) = num2str(1);
            filename(19) = num2str(1);
            
            %filename(17) = num2str(5);
            filename(20) = num2str(5);
        end

        % Read the noisy speech file
        % [input_main_noisy_speech_signal_z, ~] = readwav((filename));
        
        filename99 = filename(1:4);
        
        % use: sp01_16kHz
        filename99(5) = '_';
        filename99(6) = num2str(1);
        filename99(7) = num2str(6);
        filename99(8) = 'k';
        filename99(9) = 'H';
        filename99(10) = 'z';

        % Read the noise file
        % [true, ~] = readwav((filename99));
        
        % input_noise_noiseOnly_signal = input_main_noisy_speech_signal_z - true;
        
        %[MSE_12, MSE_notRelative_12, MedSE_final_12] = main_proposed_algorithm_project( 8, input_main_noisy_speech_signal_z, input_noise_noiseOnly_signal );

        filename
        
        % filename = 'sp01_babble_sn5';
        [s,fs] = readwav(fullfile(filename));

        z = s;
        z = [z; z; z; z; z];

        %filename = 'sp01_16kHz';
        [s2,fs2] = readwav(fullfile(filename99));

        rrrqae = s - s2;
        rrrqae = [rrrqae; rrrqae; rrrqae; rrrqae; rrrqae];
        
        %[total_noise_ps, true_total_noise_ps] = new_main2_main_proposed_algorithm_project( 8, z, rrrqae  );
        fs = 16000;
        ninc=round(0.020*fs/2);   % frame increment [fs=sample frequency]
        ovf=4;                  % overlap factor
        f=rfft(enframe(z,hanning(ovf*ninc,'periodic'),ninc),ovf*ninc,2);
        f=f.*conj(f);           % convert to power spectrum

        total_noise_ps = estnoisem(f,ninc/fs);
        
        [asgas ssss] = size(total_noise_ps);
        for pps = 1 : asgas
            total_noise_ps2(pps,:) = [total_noise_ps(pps,:) flipud(total_noise_ps(pps,2:end-1))];
        end
        total_noise_ps = total_noise_ps2;
        
        % true_total_noise_ps = true_noise_parameters_my_proposed_algorithm( filename, filename99 );
        true_total_noise_ps = true2_noise_parameters_my_proposed_algorithm( rrrqae );
        
        [MSE_3, ~] = MSE_function(true_total_noise_ps, total_noise_ps);
                
        MSE_3_total = [MSE_3_total MSE_3];
                
    end
end

MSE_3_final_total = mean(MSE_3_total);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% Create a change in power in the signal 

% We load two signals
filename = 'sp01_train_sn0';
s = wavread(fullfile(filename));

% We load two signals
filename = 'sp01_train_sn15';
y = wavread(fullfile(filename));

% We combine the two signals
g = [y' s'];
s = g;

% Plot the result
figure; set(gcf,'Color','w');
plot(0:((90116/16000)/length(s)):(90116-((90116/length(s))))/16000,s);
set(gca,'FontSize',15);
title('The SNR changes. SNR = 15 dB and 0 dB. Time domain representation of noisy speech signal.'); 
xlabel('Time in Seconds (seconds)');
ylabel('Amplitude');
grid on; grid minor;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% Implement the proposed algorithm

% Initialize the matrix
total_noise_power_frame = [];

count1 = 0;
count2 = 0;

% Define the file name
%filename = 'sp01_train_sn10';
%filename = 'sp01_train_sn0';
%filename = 'long_file_noChangeSNR';
filename = 'long_signal_changeSNR';
%filename = 'sp01_car_sn0';

% We initialize the matrix
total_speech_active_frames_noisy = [];

% Read the .wav file
[x, Srate, bits] = wavread(filename);	

% % Read the .wav file
% filename = 'sp01_babble_sn0';
% [x2, Srate, bits] = wavread(filename);	
% 
% x = [x' x2'];
% x = x';

% Define the frame size in samples
len=floor(20*Srate/1000); 
if rem(len,2)==1, len=len+1; end;

% window overlap in percent of frame size
PERC=50; 
len1=floor(len*PERC/100);
len2=len-len1;

% define window
win=hanning(len);  

% normalize window for equal level output 
win = win*len2/sum(win);  

% Noise magnitude calculations - assuming that the first 6 frames is noise/silence
nFFT=2*len;
j=1;
noise_mean=zeros(nFFT,1);
for k=1:6
    noise_mean=noise_mean+abs(fft(win.*x(j:j+len-1),nFFT));
    j=j+len;
end
noise_mu=noise_mean/6;
noise_mu2=noise_mu.^2;

%--- allocate memory and initialize various variables
k=1;
img=sqrt(-1);

x_old=zeros(len1,1);
n_old = zeros(len1,1);

Nframes=floor(length(x)/len2)-1;
xfinal=zeros(Nframes*len2,1);

% --------------- Initialize parameters ------------
k=1;
aa=0.98;
eta= 0.15;
mu=0.98;
c=sqrt(pi)/2;
qk=0.3;
qkr=(1-qk)/qk;
ksi_min=10^(-25/10); 

% Main for loop
for n=1:Nframes
    insign=win.*x(k:k+len-1);

    %--- Take fourier transform of  frame
    %
    spec=fft(insign,nFFT);
    sig=abs(spec); % compute the magnitude

    sig2=sig.^2;

    gammak=min(sig2./noise_mu2,40);  % posteriori SNR

    if n==1
        ksi=aa+(1-aa)*max(gammak-1,0);
    else
        ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     

        % decision-direct estimate of a priori SNR
        ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
    end

    log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
    vad_decision = sum(log_sigma_k)/nFFT;    

    % Noise only frame found
    if (vad_decision < eta) 
        noise_mu2 = mu * noise_mu2 + (1- mu) * sig2;
        %noise_mu2 = sig2;
    else
        % Initialize the matrix
        lamda_l_k = [];

        % for all frequency bins
        for freq_loop = 1 : length(sig2)

            % Use the likelihood ratio test
            z = 10*log10(sig2(freq_loop));
            lamda_ratio_test = (9.4970 * (exp(- (z + 67.5990)^2 / (180.7520)))) / (9.5066 * (exp(- (z + 66.8998)^2 / (180.3850))));
            %lamda_ratio_test = (9.4970 * (exp(- (z)^2 / (180.7520)))) / (9.5066 * (exp(- (z)^2 / (180.3850))));
           
            % Store the result
            lamda_l_k = [lamda_l_k lamda_ratio_test];
        end

        % Initialize the matrices
        lamda_l = 0;
        L = length(lamda_l_k);

        % for loop to compute the geometric mean
        for i = 1 : L
            lamda_l = lamda_l + log(lamda_l_k(i));
        end
        lamda_l = lamda_l / L;

        % Initialize the matrix
        power_p = [];

%        if (lamda_l >= 0)
        if (lamda_l >= log(0.46/0.54))
            for freq_loop = 1 : length(sig2)
                z = 10*log10(sig2(freq_loop));
                power_p = [power_p (9.4970 * (exp(- (z + 67.5990)^2 / (180.7520))))];
                %power_p = [power_p 10^((9.4970 * (exp(- (z)^2 / (180.7520))))/10)];
            end
            count1 = count1 + 1;
%        elseif (lamda_l < 0)
        elseif (lamda_l < log(0.46/0.54))
            for freq_loop = 1 : length(sig2)
                z = 10*log10(sig2(freq_loop));
                power_p = [power_p (9.5066 * (exp(- (z + 66.8998)^2 / (180.3850))))];
                %power_p = [power_p 10^((9.5066 * (exp(- (z)^2 / (180.3850))))/10)];
            end
            count2 = count2 + 1;
        end

        % We need to find the a posteriori probability
        power_p = power_p' .* sig2;
        power_p = power(10, (power_p / 10));        
        
        %noise_mu2 = power_p;
        noise_mu2 = mu * noise_mu2 + (1- mu) * power_p;
     end

    % Define "noise_mu"
    noise_mu = sqrt(noise_mu2);
    
    % Store the noise power
    total_noise_power_frame = [total_noise_power_frame; noise_mu2'];
    
    vk=ksi.*gammak./(1+ksi);
    j0 = besseli(0,vk/2);
    j1 = besseli(1,vk/2);

    C=exp(-0.5*vk);
    A=((c*(vk.^0.5)).*C)./gammak;
    B=(1+vk).*j0+vk.*j1;
    hw=A.*B;

    sig=sig.*hw;

    % save for estimation of a priori SNR in next frame
    Xk_prev=sig.^2;  

    xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);

    xi_w= real( xi_w);

   % --- Overlap and add ---------------
    xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
    x_old= xi_w(len1+ 1: len);

   % --- Overlap and add ---------------
   ni_time = real(ifft(noise_mu));
   ni_time_final(k:k+len2-1) = n_old + ni_time(1:len1);
   n_old = ni_time(1+len1:len);

   k=k+len2;
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Plot the noise 
[t, ~, b] = spgrambw(ni_time_final, 16000, 'd');
figure; set(gcf,'Color','w');
plot(t,b(:,80));
grid on; grid minor;
set(gca,'FontSize',15);
title('Power in dB at 500 Hz.'); 
xlabel('Time in Seconds (seconds)'); ylabel('Power (dB)');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Load the noise
%test_noise_power1 = load('noisePowerSpectrum_SNR_10_train.mat');
%test_noise_power1 = load('noisePowerSpectrum_SNR_0_train.mat');
%test_noise_power1 = load('long_signal_noisePowerSpectrum_SNR_0_train.mat');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the matrices
power_at_800Hz_noisy = [];
power_at_500Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define the frequency axis
temp = size(test_noise_power1.total_noise_power_frame, 1) / 2;
% freq = 0:Fs/temp:Fs/2;
freq = 0:8000/320:8000-(8000/320);

% Define the middle
middle_d =  size(test_noise_power1.total_noise_power_frame, 2) - (size(test_noise_power1.total_noise_power_frame, 2)/2) + 1;

for i = 1 : size(test_noise_power1.total_noise_power_frame, 1)
    % Define the PSD
    psdx = test_noise_power1.total_noise_power_frame(i, middle_d:end);
    
    % plot(freq,10*log10(psdx(end-length(freq)+1:end))); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the power at 800 Hz
    idx = find(freq==800);
    power_at_800Hz_noisy = [power_at_800Hz_noisy 10*log10(psdx(idx))];

    % Find the power at 500 Hz
    idx = find(freq==500);
    power_at_500Hz_noisy = [power_at_500Hz_noisy 10*log10(psdx(idx))];
end

% % We plot the Power Vs Time graph
% figure; set(gcf,'Color','w');
% plot(0.01 : 0.01 : length(power_at_800Hz_noisy)*0.01, power_at_800Hz_noisy)
% grid on; grid minor;
% set(gca,'FontSize',15);
% title('Power in dB at 800 Hz.'); 
% xlabel('Time in Seconds (seconds)'); ylabel('Power (dB)');

% We plot the Power Vs Time graph
figure; set(gcf,'Color','w');
plot(0.01 : 0.01 : length(power_at_500Hz_noisy)*0.01, power_at_500Hz_noisy)
grid on; grid minor;
set(gca,'FontSize',15);
title('Power in dB at 500 Hz.'); 
xlabel('Time in Seconds (seconds)'); ylabel('Power (dB)');
xlim([0.11 11.25]);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------











%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Compare with other algorithms 
% Continuous minimal tracking, Doblinger (1995), (use: 'doblinger')  

% Define the method to be used 
method = 'doblinger';

% Call the function so as to find the noise power spectrum for each frame
%[doblinger_noise_ps] = noise_parameters('long_signal_changeSNR', method);
[doblinger_noise_ps] = noise_parameters('long_file_noChangeSNR', method);

% Define the matrices
power_at_500Hz_noisy = [];

% Define Fs
Fs = 16000;

% Define the frequency axis
temp = size(doblinger_noise_ps, 1) / 2;
% freq = 0:Fs/temp:Fs/2;
freq = 0:8000/320:8000-(8000/320);

% Define the middle
middle_d =  size(doblinger_noise_ps, 2) - (size(doblinger_noise_ps, 2)/2) + 1;

for i = 1 : size(doblinger_noise_ps, 1)
    % Define the PSD
    psdx = doblinger_noise_ps(i, middle_d:end);
    
    % plot(freq,10*log10(psdx(end-length(freq)+1:end))); 
    % grid on; title('Periodogram Using FFT');
    % xlabel('Frequency (Hz)');
    % ylabel('Power/Frequency (dB/Hz)');

    % Find the power at 500 Hz
    idx = find(freq==500);
    power_at_500Hz_noisy = [power_at_500Hz_noisy 10*log10(psdx(idx))];
end

% We plot the Power Vs Time graph
figure; set(gcf,'Color','w');
plot(0.01 : 0.01 : length(power_at_500Hz_noisy)*0.01, power_at_500Hz_noisy)
grid on; grid minor;
set(gca,'FontSize',15);
title('Power in dB at 500 Hz.'); 
xlabel('Time in Seconds (seconds)'); ylabel('Power (dB)');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

plot(doblinger_noise_ps(21,:)); hold on; 
plot(test_noise_power1.total_noise_power_frame(21,:)-0.9412, 'r'); hold off;

% plot(10*log10(doblinger_noise_ps(21,:))); hold on; 
% a = test_noise_power1.total_noise_power_frame(21,:);
% plot(10*log10(a), 'r'); hold off;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------













%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%filename = 'long_file_noChangeSNR';
filename = 'long_signal_changeSNR';
figure; set(gcf,'Color','w');
h1 = specsub_ns2(filename,'martin','asdfasaaa');
% Set the color of the plot
set(h1, 'Color', [0 0 1]);
%set(h1, 'Color', [1 0 0]);

set(gca,'FontSize',15);
title('The SNR changes. Power Vs Time plot of the noisy speech signal.'); 
xlabel('Time in Seconds (seconds)');
ylabel('Power (dB)');
grid on; grid minor;
hold on;

filename = 'long_signal_changeSNR';
h2 = specsub_ns2(filename,'mcra','asdfasaaa');
% Set the color of the plot
set(h2, 'Color', [1 0 0]);

set(gca,'FontSize',15);
grid on; grid minor;
hold on;

filename = 'long_signal_changeSNR';
h2 = specsub_ns2(filename,'mcra2','asdfasaaa');
% Set the color of the plot
set(h2, 'Color', [0 1 0]);

set(gca,'FontSize',15);
grid on; grid minor;
hold off;
legend('Martin MS', 'MCRA', 'MCRA2');

% filename = 'long_signal_changeSNR';
% s = wavread(fullfile(filename));
% figure; set(gcf,'Color','w');
% plot(0:((360464/16000)/length(s)):(360464-((360464/length(s))))/16000,s);
% set(gca,'FontSize',15);
% title('The SNR changes. Time domain representation of noisy speech signal.'); 
% xlabel('Time in Seconds (seconds)');
% ylabel('Amplitude');
% grid on; grid minor;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------













%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the filename of the .wav file
%filename = 'sp01_white_noise_15dB';
filename = 'long_signal_changeSNR';
%filename = 'sp01_airport_sn0';

% Read the .wav file
[x, Srate, ~] = wavread(filename);

% Call the "noisePowProposed" function of T.Gerkmann and R.C.Hendriks
noisePowMat = noisePowProposed(x, Srate);
% noisePowMat has dimensions: (noise psd) x (frames)
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Paper: "MMSE BASED NOISE PSD TRACKING WITH LOW COMPLEXITY"
% Author: T.Gerkmann and R.C.Hendriks

% Add datapath
addpath ./TabGenGam

% Define the filename of the .wav file
%filename = 'sp01_white_noise_15dB';
filename = 'long_signal_changeSNR';
%filename = 'sp01_airport_sn0';

% Read the .wav file
[x, Srate, ~] = wavread(filename);

% Call the "noisePowProposed" function of T.Gerkmann and R.C.Hendriks
[shat, noise_psd_matrix, T] = noise_psd_tracker(x, Srate);
% noisePowMat has dimensions: (noise psd) x (frames)
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% Implement the proposed algorithm

% Call the proposed function
%filename = 'long_signal_changeSNR';
filename = 'sp01_airport_sn0';

[power_noise, x_final] = noise_my(filename, 1, 0.9, 0.6, 0.8);
%soundsc(x_final, 16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

filename = 'sp01_airport_sn15';
s = readwav(fullfile(filename));
fs = 16000;

filename = 'sp01';
[s2,fs2] = readwav(fullfile(filename));
f = resample(s2,2,1);

% Define the noise signal
noise_signal = s-f;

% Call the proposed function
filename = 'sp01_airport_sn15';
power_noise = noise_my(filename, 1, 0.8, 0.6, 0.7);
%[power_noise, x_final] = noise_my3(filename, 1, 0.9, 0.6, 0.8);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Plot the noise 
[t, ~, b] = spgrambw(noise_signal, 16000, 'd');
%figure; set(gcf,'Color','w');
hold on;
plot(t,b(:,80), 'r');
grid on; grid minor;
set(gca,'FontSize',15);
title(['SNR = 15 dB. Power in dB at 500 Hz. The mu parameter is adaptively changed.']); 
xlabel('Time in Seconds (seconds)'); ylabel('Power (dB)');
legend('Estimated Noise','True Noise');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Find the true value
total_noise_power_frame = test_proposedAlgorithm(noise_signal);
% We use "total_noise_power_frame" and "power_noise"

desired = [];
for i = 1 : size(total_noise_power_frame,1)
    desired = [desired total_noise_power_frame(i,:)];
end

observed = [];
for i = 1 : size(power_noise,1)
    observed = [observed power_noise(i,:)];
end

% MSE
MSE = mean((desired - observed).^2);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% RMS error
x_algorithm = observed;
x_reference = desired;

%RMS error: for loop
RMS_error_1 = 0;
n = length(x_reference);

%for loop
for i = 1 : n
    RMS_error_1 = RMS_error_1 + ((x_reference(i) - x_algorithm(i))^2);
end

RMS_error_1 = sqrt((1/n) * RMS_error_1);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Repeat the above procedure
RMS_total =[];
MSE_total = [];
MSE_f = 0;
Nframes_total = [];

% Repeat the above procedure
desired_total = [];
observed_total = [];

%filename = 'sp01_airport_sn0';
filename = 'sp01_exhibition_sn0';

% Outter for loop
for outter_loop = 1 : 1
    for inner_loop = 1 : 30
        % We change the filename based on inner_loop
        % We use a temporary variable
        temp_variable = num2str(inner_loop);

        % We now change the filename based on inner_loop
        if (inner_loop < 10) 
            filename(3) = num2str(0);
            filename(4) = num2str(temp_variable(1));
        else
            filename(3) = num2str(temp_variable(1));
            filename(4) = num2str(temp_variable(2));
        end

        % If statement to change the filename
        if (outter_loop == 1) 
            filename(19) = num2str(0);
        elseif (outter_loop == 2) 
            filename(19) = num2str(5);
        elseif (outter_loop == 3) 
            filename(19) = num2str(1);
            filename(20) = num2str(0);
        else 
            filename(19) = num2str(1);
            filename(20) = num2str(5);
        end
        
        if (strcmp(filename, 'sp04_babble_sn10') ~= 1)
        s = readwav(fullfile(filename));
        fs = 16000;

        if (inner_loop < 10) 
            filename2 = strcat('sp', num2str(0), num2str(inner_loop));
        else
            filename2 = strcat('sp', num2str(inner_loop));
        end
        
        s2 = readwav(fullfile(filename2));
        f = resample(s2,2,1);

        % Define the noise signal
        noise_signal = s-f;

        % Call the proposed function
        power_noise = noise_my(filename, 0, 0.7, 0.3, 0.7);
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------

        % Find the true value
        [total_noise_power_frame, Nframes] = test_proposedAlgorithm(noise_signal);
        % We use "total_noise_power_frame" and "power_noise"
        
        Nframes_total = [Nframes_total Nframes];
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------

        desired = [];
        for i = 1 : size(total_noise_power_frame,1)
            desired = [desired total_noise_power_frame(i,:)];
        end

        observed = [];
        for i = 1 : size(power_noise,1)
            observed = [observed power_noise(i,:)];
        end

        observed_total = [observed_total observed];
        desired_total = [desired_total desired];
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------

        % MSE
        MSE = mean((desired - observed).^2);
        MSE_total = [MSE_total MSE];
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------
        
        x_algorithm = observed;
        x_reference = desired;

        %RMS error: for loop
        RMS_error_1 = 0;
        n = length(x_reference);

        %for loop
        for i = 1 : n
            RMS_error_1 = RMS_error_1 + ((x_reference(i) - x_algorithm(i))^2);
        end

        RMS_error_1 = sqrt((1/n) * RMS_error_1);
        RMS_total = [RMS_total RMS_error_1];
    end
    end
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

x_reference_breaths = desired_total;
x_algorithm = observed_total;

%find the values of the 2 axes
mean_x_axis = (1/2) * (x_reference_breaths + x_algorithm);
difference_y_axis = x_reference_breaths - x_algorithm;


%plot the results
figure; set(gcf,'Color','w');
plot(mean_x_axis, difference_y_axis, 'o'); hold on;

%plot the horizontal lines
x = [0, 14*10^6];
y = [mean(difference_y_axis), mean(difference_y_axis)];
plot(x,y, 'r'); hold on;
%plot the horizontal lines
x = [0, 14*10^6];
y = [mean(difference_y_axis)+(2*std(difference_y_axis)), mean(difference_y_axis)+(2*std(difference_y_axis))];
plot(x,y, 'r'); hold on;
%plot the horizontal lines
x = [0, 14*10^6];
y = [mean(difference_y_axis)-(2*std(difference_y_axis)), mean(difference_y_axis)-(2*std(difference_y_axis))];
plot(x,y, 'r'); hold off;
grid on; grid minor; 

%name the axes
set(gca,'FontSize',15);
title(['Bland-Altman plot for all the frames of the sentences with exhibition noise at SNR = 0 dB.']);
xlabel('Mean of the 2 measurements'); 
ylabel('Difference of the 2 measurements');















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

filename = 'sp01_airport_sn0';
s = readwav(fullfile(filename));
fs = 16000;

filename = 'sp01';
[s2,fs2] = readwav(fullfile(filename));
f = resample(s2,2,1);

% Define the noise signal
noise_signal = s-f;

% Find the true value
total_noise_power_frame = test_proposedAlgorithm(noise_signal);
% We use "total_noise_power_frame" and "power_noise"

% % Call the proposed function
% filename = 'sp01_airport_sn0';
% [~, ~, MSE] = noise_my2(filename, 0, 0.7, 0.3, 0.7, total_noise_power_frame);

% Call the "noisePowProposed" function of T.Gerkmann and R.C.Hendriks
filename = 'sp01_airport_sn0';
[x, Srate, ~] = wavread(filename);
[~, MSE2] = noisePowProposed2(x, Srate, total_noise_power_frame);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------







%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

filename = 'sp01_airport_sn0';
[power_noise, x_final, counter_pp] = noise_my3(filename, 0, 0.7, 0.3, 0.7);
%soundsc(x_final, 16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------












%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Create a GUI table for R.Martin's MS algorithm

% We use R.Martin's MS algorithm
method = 'martin';

filename = 'sp04_babble_sn10';
outfile = strcat('output_', filename);

noise_ps_martin = specsub_ns(filename, method, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp04';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = pesq_val(2);

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp04_babble_sn5';

outfile = strcat('output_', filename);

%method = 'mcra';
method = 'martin';

noise_ps_martin = specsub_ns(filename, method, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp04';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = pesq_val(2);

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final2 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp04_babble_sn15';

outfile = strcat('output_', filename);

method = 'martin';

noise_ps_martin = specsub_ns(filename, method, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp04';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = pesq_val(2);

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final3 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp04_babble_sn0';

outfile = strcat('output_', filename);

method = 'martin';

noise_ps_martin = specsub_ns(filename, method, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp04';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = pesq_val(2);

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final4 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

%create a table with the values calculated
T = table(final4, final2, final3, final);
T.Properties.VariableNames = {'TestResult', 'TestResult2', 'TestResult3', 'TestResult4'};

%print the results
disp(['The following table shows the results that we just calculated.']);
disp(T);
disp(['   ']); disp(['   ']);

%set the figure data
f = figure; set(gcf,'Color','w');
set(f, 'Position', [500 400 300 100]);
dat =  {T{1,1}, T{1,2}, T{1,3}, T{1,4};...
        T{2,1}, T{2,2}, T{2,3}, T{2,4};...   
        T{3,1}, T{3,2}, T{3,3}, T{3,4};...
        T{4,1}, T{4,2}, T{4,3}, T{4,4};...
        T{5,1}, T{5,2}, T{5,3}, T{5,4};...
        T{6,1}, T{6,2}, T{6,3}, T{6,4};...
        T{7,1}, T{7,2}, T{7,3}, T{7,4};...
        T{8,1}, T{8,2}, T{8,3}, T{8,4};...
        T{9,1}, T{9,2}, T{9,3}, T{9,4};...
        T{10,1}, T{10,2}, T{10,3}, T{10,4};...
        T{11,1}, T{11,2}, T{11,3}, T{11,4};...
        T{12,1}, T{12,2}, T{12,3}, T{12,4};...
        T{13,1}, T{13,2}, T{13,3}, T{13,4};...
        T{14,1}, T{14,2}, T{14,3}, T{14,4};...
        T{15,1}, T{15,2}, T{15,3}, T{15,4};...
        T{16,1}, T{16,2}, T{16,3}, T{16,4};...
        T{17,1}, T{17,2}, T{17,3}, T{17,4};...
        T{18,1}, T{18,2}, T{18,3}, T{18,4};};

%set the table as a figure
columnname =   {'Results for 15 dB SNR', 'Results for 10 dB SNR', 'Results for 5 dB SNR', 'Results for 0 dB SNR'};
rnames = {'SNR mean', 'Segmental SNR', 'WSS mean', 'LLR mean', 'IS mean', 'CEP mean', ...
    'fw SNR segmental', 'SIG' 'BAK', 'OVL', 'SIG2' 'BAK2', 'OVL2', 'raw PESQ', 'MOS-valued PESQ', 'Csig', 'Cbak', 'Covl'};
columnformat = {'numeric'}; 
t = uitable('Units','normalized','Position',...
            [0.05 0.05 0.755 0.87], 'Data', dat,... 
            'ColumnName', columnname,...
            'RowName', rnames, ...
            'ColumnFormat', columnformat); 
        
%set the title
%set(gca,'FontSize',15); 
title('Table: R.Martin MS Algorithm. Only babble noise is used. Only one speaker is tested.');        
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Create a GUI table for the proposed algorithm

% We use the proposed algorithm
filename = 'sp01_babble_sn10';
outfile = strcat('output_', filename);

%noise_ps_martin = specsub_ns(filename, method, outfile);
[noise_ps_martin, xfinal] = noise_my(filename, 0, 0.7, 0.3, 0.7);
writewav(xfinal, 16000, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp01_16kHz';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = 0;
if (length(pesq_val) > 1)
    MOS_mapped_PESQ = pesq_val(2);
end

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp01_babble_sn5';
outfile = strcat('output_', filename);

%noise_ps_martin = specsub_ns(filename, method, outfile);
[noise_ps_martin, xfinal] = noise_my(filename, 0, 0.7, 0.3, 0.7);
wavwrite(xfinal, 16000, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp01_16kHz';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = 0;
if (length(pesq_val) > 1)
    MOS_mapped_PESQ = pesq_val(2);
end

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final2 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp01_babble_sn15';
outfile = strcat('output_', filename);

%noise_ps_martin = specsub_ns(filename, method, outfile);
[noise_ps_martin, xfinal] = noise_my(filename, 0, 0.7, 0.3, 0.7);
wavwrite(xfinal, 16000, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp01_16kHz';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = 0;
if (length(pesq_val) > 1)
    MOS_mapped_PESQ = pesq_val(2);
end

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final3 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp01_babble_sn0';
outfile = strcat('output_', filename);

%noise_ps_martin = specsub_ns(filename, method, outfile);
[noise_ps_martin, xfinal] = noise_my(filename, 0, 0.7, 0.3, 0.7);
wavwrite(xfinal, 16000, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp01_16kHz';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = 0;
if (length(pesq_val) > 1)
    MOS_mapped_PESQ = pesq_val(2);
end

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final4 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

%create a table with the values calculated
T = table(final4, final2, final3, final);
T.Properties.VariableNames = {'TestResult', 'TestResult2', 'TestResult3', 'TestResult4'};

%print the results
disp(['The following table shows the results that we just calculated.']);
disp(T);
disp(['   ']); disp(['   ']);

%set the figure data
f = figure; set(gcf,'Color','w');
set(f, 'Position', [500 400 300 100]);
dat =  {T{1,1}, T{1,2}, T{1,3}, T{1,4};...
        T{2,1}, T{2,2}, T{2,3}, T{2,4};...   
        T{3,1}, T{3,2}, T{3,3}, T{3,4};...
        T{4,1}, T{4,2}, T{4,3}, T{4,4};...
        T{5,1}, T{5,2}, T{5,3}, T{5,4};...
        T{6,1}, T{6,2}, T{6,3}, T{6,4};...
        T{7,1}, T{7,2}, T{7,3}, T{7,4};...
        T{8,1}, T{8,2}, T{8,3}, T{8,4};...
        T{9,1}, T{9,2}, T{9,3}, T{9,4};...
        T{10,1}, T{10,2}, T{10,3}, T{10,4};...
        T{11,1}, T{11,2}, T{11,3}, T{11,4};...
        T{12,1}, T{12,2}, T{12,3}, T{12,4};...
        T{13,1}, T{13,2}, T{13,3}, T{13,4};...
        T{14,1}, T{14,2}, T{14,3}, T{14,4};...
        T{15,1}, T{15,2}, T{15,3}, T{15,4};...
        T{16,1}, T{16,2}, T{16,3}, T{16,4};...
        T{17,1}, T{17,2}, T{17,3}, T{17,4};...
        T{18,1}, T{18,2}, T{18,3}, T{18,4};};

%set the table as a figure
columnname =   {'Results for 15 dB SNR', 'Results for 10 dB SNR', 'Results for 5 dB SNR', 'Results for 0 dB SNR'};
rnames = {'SNR mean', 'Segmental SNR', 'WSS mean', 'LLR mean', 'IS mean', 'CEP mean', ...
    'fw SNR segmental', 'SIG' 'BAK', 'OVL', 'SIG2' 'BAK2', 'OVL2', 'raw PESQ', 'MOS-valued PESQ', 'Csig', 'Cbak', 'Covl'};
columnformat = {'numeric'}; 
t = uitable('Units','normalized','Position',...
            [0.05 0.05 0.755 0.87], 'Data', dat,... 
            'ColumnName', columnname,...
            'RowName', rnames, ...
            'ColumnFormat', columnformat); 
        
%set the title
%set(gca,'FontSize',15); 
title('Table: Proposed Algorithm. Only babble noise is used. Only one speaker is tested.');        
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use Spectral Subtraction (SS)
xfinal = noise_my_SS('sp01_restaurant_sn0', 0.7, 0.3, 0.7);
%soundsc(xfianl, 16000);

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Create a GUI table for the proposed algorithm

% We use the proposed algorithm
filename = 'sp01_babble_sn10';
outfile = strcat('output_', filename);

%noise_ps_martin = specsub_ns(filename, method, outfile);
xfinal = noise_my_SS(filename, 0.7, 0.3, 0.7);
writewav(xfinal, 16000, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp01_16kHz';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = 0;
if (length(pesq_val) > 1)
    MOS_mapped_PESQ = pesq_val(2);
end

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp01_babble_sn5';
outfile = strcat('output_', filename);

%noise_ps_martin = specsub_ns(filename, method, outfile);
xfinal = noise_my_SS(filename, 0.7, 0.3, 0.7);
wavwrite(xfinal, 16000, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp01_16kHz';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = 0;
if (length(pesq_val) > 1)
    MOS_mapped_PESQ = pesq_val(2);
end

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final2 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp01_babble_sn15';
outfile = strcat('output_', filename);

%noise_ps_martin = specsub_ns(filename, method, outfile);
xfinal = noise_my_SS(filename, 0.7, 0.3, 0.7);
wavwrite(xfinal, 16000, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp01_16kHz';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = 0;
if (length(pesq_val) > 1)
    MOS_mapped_PESQ = pesq_val(2);
end

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final3 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

filename = 'sp01_babble_sn0';
outfile = strcat('output_', filename);

%noise_ps_martin = specsub_ns(filename, method, outfile);
xfinal = noise_my_SS(filename, 0.7, 0.3, 0.7);
wavwrite(xfinal, 16000, outfile);

[s, fs] = wavread(outfile);
%soundsc(s, fs);

filename2 = 'sp01_16kHz';

[snr_mean, segsnr_mean]= comp_snr(filename2, outfile);
%      where 'snr_mean' is the global overall SNR and 'segsnr_mean' is the 
%      segmental SNR.

wss_mean = comp_wss(filename2, outfile);

llr_mean= comp_llr(filename2, outfile);

is_mean = comp_is(filename2, outfile);

cep_mean = comp_cep(filename2, outfile);

fwSNRseg = comp_fwseg(filename2, outfile);

[SIG,BAK,OVL] = comp_fwseg_variant(filename2, outfile);
%	where   'SIG' is the predicted rating of speech distortion,
%		'BAK' is the predicted rating of background noise distortion,
%		'OVL' is the predicted rating of overall quality.

[SIG2,BAK2,OVL2] = comp_fwseg_mars(filename2, outfile);

pesq_val = pesq(strcat(filename2, '.wav'), strcat(outfile, '.wav')); 
raw_PESQ = pesq_val(1);
MOS_mapped_PESQ = 0;
if (length(pesq_val) > 1)
    MOS_mapped_PESQ = pesq_val(2);
end

[Csig,Cbak,Covl] = composite(strcat(filename2, '.wav'), strcat(outfile, '.wav'));
%	where   'Csig' is the predicted rating of speech distortion,
%		'Cbak' is the predicted rating of background noise distortion,
%		'Covl' is the predicted rating of overall quality.

final4 = final;

%test = ['SNRmean', 'SegSNR', 'WSSmean', 'LLRmean', 'ISmean', 'CEPmean', ...
%    'fwSNRseg', 'SIG' 'BAK', 'OVL', 'SIG' 'BAK', 'OVL', 'PESQ', 'Csig', 'Cbak', 'Covl'];
final = [snr_mean, segsnr_mean,wss_mean,llr_mean,is_mean,cep_mean,fwSNRseg,SIG,BAK,OVL,SIG2,BAK2,OVL2,raw_PESQ,MOS_mapped_PESQ,Csig,Cbak,Covl];
final = final';

%create a table with the values calculated
T = table(final4, final2, final3, final);
T.Properties.VariableNames = {'TestResult', 'TestResult2', 'TestResult3', 'TestResult4'};

%print the results
disp(['The following table shows the results that we just calculated.']);
disp(T);
disp(['   ']); disp(['   ']);

%set the figure data
f = figure; set(gcf,'Color','w');
set(f, 'Position', [500 400 300 100]);
dat =  {T{1,1}, T{1,2}, T{1,3}, T{1,4};...
        T{2,1}, T{2,2}, T{2,3}, T{2,4};...   
        T{3,1}, T{3,2}, T{3,3}, T{3,4};...
        T{4,1}, T{4,2}, T{4,3}, T{4,4};...
        T{5,1}, T{5,2}, T{5,3}, T{5,4};...
        T{6,1}, T{6,2}, T{6,3}, T{6,4};...
        T{7,1}, T{7,2}, T{7,3}, T{7,4};...
        T{8,1}, T{8,2}, T{8,3}, T{8,4};...
        T{9,1}, T{9,2}, T{9,3}, T{9,4};...
        T{10,1}, T{10,2}, T{10,3}, T{10,4};...
        T{11,1}, T{11,2}, T{11,3}, T{11,4};...
        T{12,1}, T{12,2}, T{12,3}, T{12,4};...
        T{13,1}, T{13,2}, T{13,3}, T{13,4};...
        T{14,1}, T{14,2}, T{14,3}, T{14,4};...
        T{15,1}, T{15,2}, T{15,3}, T{15,4};...
        T{16,1}, T{16,2}, T{16,3}, T{16,4};...
        T{17,1}, T{17,2}, T{17,3}, T{17,4};...
        T{18,1}, T{18,2}, T{18,3}, T{18,4};};

%set the table as a figure
columnname =   {'Results for 15 dB SNR', 'Results for 10 dB SNR', 'Results for 5 dB SNR', 'Results for 0 dB SNR'};
rnames = {'SNR mean', 'Segmental SNR', 'WSS mean', 'LLR mean', 'IS mean', 'CEP mean', ...
    'fw SNR segmental', 'SIG' 'BAK', 'OVL', 'SIG2' 'BAK2', 'OVL2', 'raw PESQ', 'MOS-valued PESQ', 'Csig', 'Cbak', 'Covl'};
columnformat = {'numeric'}; 
t = uitable('Units','normalized','Position',...
            [0.05 0.05 0.755 0.87], 'Data', dat,... 
            'ColumnName', columnname,...
            'RowName', rnames, ...
            'ColumnFormat', columnformat); 
        
%set the title
%set(gca,'FontSize',15); 
title('Table: Proposed Algorithm. Only babble noise is used. Only one speaker is tested.');        
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------



















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Read the input file
filename2 = 'SA1.WAV'; fs = 16000;
s = readsph(fullfile(filename2));

% Read the noise file
filename = 'destroyerengine'; 
[s2, fs2] = readwav(fullfile(filename));

% Define the SNR
snr = 15;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'dopEk',s2,fs2); 

% We use this for noise only: "z2(:,2)"
z2 = v_addnoise(s,fs,snr,'dxopEk',s2,fs2);
% We use "soundsc(z2(:,2),fs)" for noise only
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Find the true power of noise
noise_signal = z2(:,2);
true_total_noise_power_frame = test_proposedAlgorithm(noise_signal);
true_total_noise_power_frame = true_total_noise_power_frame(8:end, :);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Implement the proposed algorithm
% The "mu" parameter defines how much does the noise estimate depend on the
% previous frame noise estimate. It is adaptively changed. 
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Initialize the matrix
total_noise_power_frame = [];

% We initialize the matrix
total_speech_active_frames_noisy = [];
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PEFAC algorithm
x = z; Srate = 16000;
[fx,tx,pv,fv,prv,pru] = fxpefac(x, Srate);

% find the probability
pv = pv*10^6 / max(pv*10^6);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

n = 320;
f = enframe(x,hamming(n,'periodic'),n/2);  
f = f(8:end, :);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the frame size in samples
Srate = 16000;
len = floor(20 * Srate/1000); 
if rem(len,2)==1, len=len+1; end;

% window overlap in percent of frame size
PERC = 50; 
len1 = floor(len * PERC/100);
len2 = len - len1;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;
nFFT = 2 * len;

x_old = zeros(len1,1);
%n_old = zeros(len1,1);
xfinal = zeros(Nframes*len2, 1);

k = 1;
img = sqrt(-1);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

mu1 = 0.9;
mu2 = 0.6;
mu3 = 0.7;

threshold1 = 0.2;
%threshold1 = 0.3;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the clean speech
m1 = -18.4705;
m2 = -21.4016;
v1 = 18.0409;
v2 = 49.9072;
w1 = 0.7772;
w2 = 0.2228;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the high noise power
m3 = -46.2442;
m4 = -58.1222;
v3 = 32.9259;
v4 = 70.9205;
w3 = 0.7777;
w4 = 0.2223;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the low noise power
m5 = -53.6041;
m6 = -45.0366;
v5 = 72.7698;
v6 = 27.0183;
w5 = 0.3750;
w6 = 0.6250;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

count1 = 0;
count2 = 0;
count3 = 0;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    insign = f(n,:)';

    % Take Fourier transform of  frame
    spec = fft(insign,nFFT);

    % Compute the magnitude spectrum
    sig = abs(spec); 
   
    % Define the power of the frame
    sig2 = sig .^ 2;
    
    % Assume the first frame is noise
    if (n == 1)
        noise_mu2 = sig2;
        final = noise_mu2';
    end
    
    % Silence/Unvoiced frame found
    if (pv(n) < threshold1) 
        mu = mu1;
        noise_mu2 = mu * noise_mu2 + (1 - mu) * sig2;
        %noise_mu2 = sig2;
    % Voiced frame found
    else
% count3 = count3 + 1;
% % Initialize the matrix
% lamda_l_k = [];
% 
% % for all frequency bins
% for freq_loop = 1 : length(sig2)
% 
%     % Use the likelihood ratio test
%     z = 10*log10(sig2(freq_loop));
%             pdfCompute = pdf_compute(z, m1, m2, m3, m4, v1, v2, v3, v4, w1, w2, w3, w4);
%             pdfCompute2 = pdf_compute(z, m1, m2, m5, m6, v1, v2, v5, v6, w1, w2, w5, w6);
%     %lamda_ratio_test = (9.4970 * (exp(- (z + 67.5990)^2 / (180.7520)))) / (9.5066 * (exp(- (z + 66.8998)^2 / (180.3850))));
%     %lamda_ratio_test = (9.4970 * (exp(- (z)^2 / (180.7520)))) / (9.5066 * (exp(- (z)^2 / (180.3850))));
%     lamda_ratio_test = pdfCompute2 / pdfCompute;
%     
%     % Store the result
%     lamda_l_k = [lamda_l_k lamda_ratio_test];
% end
% 
% % Initialize the matrices
% lamda_l = 0;
% L = length(lamda_l_k);
% 
% % for loop to compute the geometric mean
% for i = 1 : L
%     lamda_l = lamda_l + log(lamda_l_k(i));
% end
% lamda_l = lamda_l / L;
% 
% % Initialize the matrix
% power_p = [];
% 
% %        if (lamda_l >= 0)
% if (lamda_l >= log(0.46/0.54))
%     for freq_loop = 1 : length(sig2)
%         z = 10*log10(sig2(freq_loop));
%         power_p = [power_p (9.4970 * (exp(- (z + 67.5990)^2 / (180.7520))))];
%         %power_p = [power_p 10^((9.4970 * (exp(- (z)^2 / (180.7520))))/10)];
%     end
%     count1 = count1 + 1;
% %        elseif (lamda_l < 0)
% else
%     for freq_loop = 1 : length(sig2)
%         z = 10*log10(sig2(freq_loop));
%         power_p = [power_p (9.5066 * (exp(- (z + 66.8998)^2 / (180.3850))))];
%         %power_p = [power_p 10^((9.5066 * (exp(- (z)^2 / (180.3850))))/10)];
%     end
%     count2 = count2 + 1;
% end
        
        % Initialize the matrices
        power_p = [];
        power_p2 = [];

        % Hypothesis H_0: Transient Noise exists
        for freq_loop = 1 : length(sig2)
            z = 10*log10(sig2(freq_loop));

            %power_p = [power_p (9.4970 * (exp(- (z + 67.5990)^2 / (180.7520))))];
            pdfCompute = pdf_compute(z, m1, m2, m3, m4, v1, v2, v3, v4, w1, w2, w3, w4);
            power_p = [power_p pdfCompute];
        end

        % Hypothesis H_1: Background Noise exists
        for freq_loop = 1 : length(sig2)
            z = 10*log10(sig2(freq_loop));
            
            %power_p2 = [power_p2 (9.5066 * (exp(- (z + 66.8998)^2 / (180.3850))))];
            pdfCompute2 = pdf_compute(z, m1, m2, m5, m6, v1, v2, v5, v6, w1, w2, w5, w6);
            power_p2 = [power_p2 pdfCompute2];
        end
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------
        
        % We need to find the a posteriori probabilities
        % Hypothesis H_0: Transient Noise exists
        power_p = (power_p * 0.4) ./ (0.4 * power_p + 0.6 * power_p2); 

        % Hypothesis H_1: Background Noise exists
        power_p2 = (power_p2 * 0.6) ./ (0.4 * power_p + 0.6 * power_p2); 
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------

        % We use a discrepancy function
        % Define the variables 
        final_function = 0;
        observed_noisy_speech = sqrt(sig2);
        previous_noise_estimate = noise_mu;

        % Implement the discrepancy function
        final_function = sum(abs(observed_noisy_speech-previous_noise_estimate)) / sum(previous_noise_estimate);
        final_function = abs(final_function);
        discrepancy_measure = min(final_function, 1);
         
        if (discrepancy_measure > 0.6)
            % We increase "power_p" by 10%
            power_p = power_p * 1.1;
            
            % We decrease "power_p2" by 10%
            power_p2 = power_p2 * 0.9;
        else 
            % We increase "power_p2" by 10%
            power_p2 = power_p2 * 1.1;
            
            % We decrease "power_p" by 10%
            power_p = power_p * 0.9;
        end
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------
        
        % find final 
        %final = power_p .* (10*log10(sig2))' + power_p2 .* (10*log10(final));
        final = (power_p/(power_p+power_p2)) .* (10*log10(sig2))' + (power_p2/(power_p+power_p2)) .* (10*log10(final));
        final = power(10, (final / 10));        
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------
        
        % Find the final "noise_mu2"
        %noise_mu2 = final';
        
        if (discrepancy_measure > 0.6)
            % Update the noise power spectrum
            mu = mu2;
            noise_mu2 = mu * noise_mu2 + (1 - mu) * final';
        else 
            % Update the noise power spectrum
            mu = mu3;
            noise_mu2 = mu * noise_mu2 + (1 - mu) * final';
        end
    end
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------

    % Define "noise_mu"
    noise_mu = sqrt(noise_mu2);

    % Store the noise power
    total_noise_power_frame = [total_noise_power_frame; noise_mu2'];
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------

    % Do overlap-add 
    xi_w = real(ifft(sig .* exp(img*angle(spec)), nFFT));
    xfinal(k:k+ len2-1)= x_old + xi_w(1:len1);
    x_old= xi_w(len1+ 1: len);

%     % Do overlap-add 
%     ni_time = real(ifft(noise_mu .* exp(img*angle(spec)), nFFT));
%     ni_time_final(k:k+len2-1) = n_old + ni_time(1:len1);
%     n_old = ni_time(1+len1:len);

    % Update k
    k = k + len2;
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use MSE
MSE = 0;

% Check the result of the algorittm
% We use: "total_noise_power_frame"
for i = 1 : size(total_noise_power_frame, 1)
    for k = 1: size(total_noise_power_frame, 2)
        % Compute the MSE
        MSE = MSE + ((total_noise_power_frame(i,k) - true_total_noise_power_frame(i,k))^2 / true_total_noise_power_frame(i,k));
    end
end

% Find the final value of the MSE
MSE = MSE / (size(total_noise_power_frame, 1)*size(total_noise_power_frame, 2));
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------













%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

filename = 'sp01_airport_sn10';
s = readwav(fullfile(filename));
fs = 16000;

filename = 'sp01';
[s2,fs2] = readwav(fullfile(filename));
f = resample(s2,2,1);

% Define the noise signal
noise_signal = s-f;

% Call the proposed function
filename = 'sp01_airport_sn10';
power_noise = noise_my7(filename, 1, 0.82, 0.6, 0.7);
%[power_noise, x_final] = noise_my3(filename, 1, 0.9, 0.6, 0.8);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Plot the noise 
[t, ~, b] = spgrambw(noise_signal, 16000, 'd');
%figure; set(gcf,'Color','w');
hold on;
plot(t,b(:,80), 'r');
grid on; grid minor;

% Set the title
set(gca,'FontSize',15);
title(['SNR = 10 dB. Power in dB at 500 Hz. The mu parameter is adaptively changed.']); 
xlabel('Time in Seconds (seconds)'); ylabel('Power (dB)');
legend('Estimated Noise','True Noise');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------



















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
% Read the input file
filename2 = 'SA1.WAV'; fs = 16000;
s = readsph(fullfile(filename2));

% Read the noise file
filename = 'destroyerengine'; 
[s2, fs2] = readwav(fullfile(filename));

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'dopEk',s2,fs2); 

% We use this for noise only: "z2(:,2)"
z2 = v_addnoise(s,fs,snr,'dxopEk',s2,fs2);
% We use "soundsc(z2(:,2),fs)" for noise only
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % Read the input file
% filename2 = 'SA1.WAV'; fs = 16000;
% s = readsph(fullfile(filename2));
% 
% % Read the noise file
% filename = 'destroyerengine'; 
% [s2, fs2] = readwav(fullfile(filename));
% 
% % Define the SNR
% snr = 0;
% 
% % Add noise from column vector n() at sample frequency fn with random start
% % Sample and wrapping around as required with a vorbis cross-fade
% z = [z; v_addnoise(s,fs,snr,'dopEk',s2,fs2)]; 
% 
% % We use this for noise only: "z2(:,2)"
% z2 = [z2; v_addnoise(s,fs,snr,'dxopEk',s2,fs2)];
% % We use "soundsc(z2(:,2),fs)" for noise only
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Find the true power of noise
noise_signal = z2(:,2);
true_total_noise_power_frame = test_proposedAlgorithm(noise_signal);
true_total_noise_power_frame = true_total_noise_power_frame(8:end, :);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Implement the proposed algorithm
% The "mu" parameter defines how much does the noise estimate depend on the
% previous frame noise estimate. It is adaptively changed. 
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Initialize the matrix
total_noise_power_frame = [];

% We initialize the matrix
total_speech_active_frames_noisy = [];
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PEFAC algorithm
x = z; Srate = 16000;
[fx,tx,pv,fv] = fxpefac(x, Srate);

% We have implemented the PEFAC algorithm
% We will use the probability of being voiced: "pv"
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

n = 320;
f = enframe(x,hamming(n,'periodic'),n/2);  
f = f(8:end, :);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the frame size in samples
Srate = 16000;
len = floor(20 * Srate/1000); 
if rem(len,2)==1, len=len+1; end;

% window overlap in percent of frame size
PERC = 50; 
len1 = floor(len * PERC/100);
len2 = len - len1;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;
nFFT = 2 * len;

x_old = zeros(len1,1);
%n_old = zeros(len1,1);
xfinal = zeros(Nframes*len2, 1);

k = 1;
img = sqrt(-1);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

mu1 = 0.8;
mu2 = 0.6;
mu3 = 0.7;

threshold1 = 0.2;
%threshold1 = 0.3;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the clean speech
m1 = -18.4705;
m2 = -21.4016;
v1 = 18.0409;
v2 = 49.9072;
w1 = 0.7772;
w2 = 0.2228;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the high noise power
m3 = -46.2442;
m4 = -58.1222;
v3 = 32.9259;
v4 = 70.9205;
w3 = 0.7777;
w4 = 0.2223;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the low noise power
m5 = -53.6041;
m6 = -45.0366;
v5 = 72.7698;
v6 = 27.0183;
w5 = 0.3750;
w6 = 0.6250;
w6 = 0.6250;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

count1 = 0;
count2 = 0;
count3 = 0;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    insign = f(n,:)';

    % Take Fourier transform of  frame
    spec = fft(insign, nFFT);

    % Compute the magnitude
    sig = abs(spec); 

    % Find the power spectrum
    %sig2 = spec .* conj(spec);
    sig2 = sig .^ 2;
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------
    
    % Assume the first frame is noise
    if (n == 1)
        noise_mu2 = sig2;
        final = noise_mu2;
    end
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------
    
%     % Compute the a posteriori SNR
%     gammak = min(sig2./noise_mu2, 40);  
% 
%     if n==1
%         ksi=aa+(1-aa)*max(gammak-1,0);
%     else
%         ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
% 
%         % decision-direct estimate of a priori SNR
%         ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%     end
% 
%     log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%     vad_decision = sum(log_sigma_k)/nFFT;    
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------

    % Noise only frame found
    %if (vad_decision < eta)
    if (pv(n) < 0.5)
        %mu = 0.9;
        mu = mu1;
        %noise_mu2 = mu * noise_mu2 + (1 - mu) * sig2;
        noise_mu2 = mu * sig2;
    else
        % Initialize the matrices
        power_p = [];
        power_p2 = [];

       % Hypothesis H_0: Transient Noise exists
       for freq_loop = 1 : length(sig2)
            z = 10*log10(sig2(freq_loop));

            %power_p = [power_p (9.4970 * (exp(- (z + 67.5990)^2 / (180.7520))))];
            pdfCompute = pdf_compute(z, m1, m2, m3, m4, v1, v2, v3, v4, w1, w2, w3, w4);
            power_p = [power_p pdfCompute];
        end

        % Hypothesis H_1: Background Noise exists
        for freq_loop = 1 : length(sig2)
            z = 10*log10(sig2(freq_loop));

            %power_p2 = [power_p2 (9.5066 * (exp(- (z + 66.8998)^2 / (180.3850))))];
            pdfCompute2 = pdf_compute(z, m1, m2, m5, m6, v1, v2, v5, v6, w1, w2, w5, w6);
            power_p2 = [power_p2 pdfCompute2];
        end
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------

        % We use a discrepancy function
        % Define the variables 
        final_function = 0;
        observed_noisy_speech = sqrt(sig2);
        previous_noise_estimate = noise_mu;

        % Implement the discrepancy function
        final_function = sum(abs(observed_noisy_speech-previous_noise_estimate)) / sum(previous_noise_estimate);
        final_function = abs(final_function);
        discrepancy_measure = min(final_function, 1);

        if (discrepancy_measure > 0.8)
            % We increase "power_p" by 5%
            power_p = power_p * 1.05;

            % We decrease "power_p2" by 5%
            power_p2 = power_p2 * 0.95;
        else 
            % We increase "power_p2" by 5%
            power_p2 = power_p2 * 1.05;

            % We decrease "power_p" by 5%
            power_p = power_p * 0.95;
        end
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------

        lamda_l_k_frameBin = power_p2 ./ power_p;

        % Use the geometric mean
        % find lamda_l_frame
        lamda_l_frame = 0;
        for u = 1 : length(sig2)
           lamda_l_frame = lamda_l_frame * lamda_l_k_frameBin(u); 
        end

        % find lamda_l_frame
        lamda_l_frame = lamda_l_frame^(1 / length(sig2));
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------

        r_aPriori = 0.6 / 0.4; 

        % find final 
        final = ((r_aPriori * lamda_l_frame) / (1 + (r_aPriori * lamda_l_frame))) * (10*log10(final))' + (1 / (1 + (r_aPriori * lamda_l_frame))) * (10*log10(sig2))'; 
        final = final';
        final = power(10, (final / 10));        
        %-------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------------------------------------------------------

        % We use the discrepancy function
        if (discrepancy_measure > 0.8)
            % Update the noise power spectrum
            %mu = 0.2;
            mu = mu2;
            %noise_mu2 = mu * noise_mu2 + (1 - mu) * final;
            noise_mu2 = (1-mu) * final;
        else
            % Update the noise power spectrum
            %mu = 0.9;
            mu = mu3;
            %noise_mu2 = mu * noise_mu2 + (1 - mu) * final;
            noise_mu2 = (1-mu) * final;
        end
    end
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------

    % Define "noise_mu"
    noise_mu = sqrt(noise_mu2);

    % Store the noise power
    total_noise_power_frame = [total_noise_power_frame; noise_mu2'];
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------

%     vk=ksi.*gammak./(1+ksi);
%     j0 = besseli(0,vk/2);
%     j1 = besseli(1,vk/2);
% 
%     C=exp(-0.5*vk);
%     A=((c*(vk.^0.5)).*C)./gammak;
%     B=(1+vk).*j0+vk.*j1;
%     hw=A.*B;
% 
%     sig=sig.*hw;
% 
%     % save for estimation of a priori SNR in next frame
%     Xk_prev=sig.^2;  
%     %-------------------------------------------------------------------------------------------------------------------------------
%     %-------------------------------------------------------------------------------------------------------------------------------
% 
%     xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%     xi_w= real( xi_w);
%     %-------------------------------------------------------------------------------------------------------------------------------
%     %-------------------------------------------------------------------------------------------------------------------------------
% 
%     % --- Overlap and add ---------------
%     xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%     x_old= xi_w(len1+ 1: len);
% 
%     % --- Overlap and add ---------------
%     ni_time = real(ifft(noise_mu));
%     ni_time_final(k:k+len2-1) = n_old + ni_time(1:len1);
%     n_old = ni_time(1+len1:len);
% 
%     % Update k
%     k = k + len2;
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use MSE
MSE = 0;

% Check the result of the algorittm
% We use: "total_noise_power_frame"
for i = 1 : size(total_noise_power_frame, 1)
    for k = 1: size(total_noise_power_frame, 2)
        % Compute the MSE
        MSE = MSE + ((total_noise_power_frame(i,k) - true_total_noise_power_frame(i,k))^2 / true_total_noise_power_frame(i,k));
    end
end

% Find the final value of the MSE
MSE = MSE / (size(total_noise_power_frame, 1)*size(total_noise_power_frame, 2));
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
% Consider the use of sub-frames

% Define the number of sub-frames
number_subFrames = 4;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Read the input file
filename2 = 'SA1.WAV'; fs = 16000;
s = readsph(fullfile(filename2));

% Read the noise file
filename = 'destroyerengine'; 
[s2, fs2] = readwav(fullfile(filename));

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'dopEk',s2,fs2); 

% We use this for noise only: "z2(:,2)"
z2 = v_addnoise(s,fs,snr,'dxopEk',s2,fs2);
% We use "soundsc(z2(:,2),fs)" for noise only
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % Read the input file
% filename2 = 'SA1.WAV'; fs = 16000;
% s = readsph(fullfile(filename2));
% 
% % Read the noise file
% filename = 'destroyerengine'; 
% [s2, fs2] = readwav(fullfile(filename));
% 
% % Define the SNR
% snr = 0;
% 
% % Add noise from column vector n() at sample frequency fn with random start
% % Sample and wrapping around as required with a vorbis cross-fade
% z = [z; v_addnoise(s,fs,snr,'dopEk',s2,fs2)]; 
% 
% % We use this for noise only: "z2(:,2)"
% z2 = [z2; v_addnoise(s,fs,snr,'dxopEk',s2,fs2)];
% % We use "soundsc(z2(:,2),fs)" for noise only
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Find the true power of noise
noise_signal = z2(:,2);
true_total_noise_power_frame = test_proposedAlgorithm(noise_signal);
true_total_noise_power_frame = true_total_noise_power_frame(8:end, :);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Implement the proposed algorithm
% The "mu" parameter defines how much does the noise estimate depend on the
% previous frame noise estimate. It is adaptively changed. 
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Initialize the matrix
total_noise_power_frame = [];

% We initialize the matrix
total_speech_active_frames_noisy = [];
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PEFAC algorithm
x = z; Srate = 16000;
[fx,tx,pv,fv] = fxpefac(x, Srate);

% We have implemented the PEFAC algorithm
% We will use the probability of being voiced: "pv"
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

n = 320;
f = enframe(x,hamming(n,'periodic'),n/2);  
f = f(8:end, :);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the frame size in samples
Srate = 16000;
len = floor(20 * Srate/1000); 
if rem(len,2)==1, len=len+1; end;

% window overlap in percent of frame size
PERC = 50; 
len1 = floor(len * PERC/100);
len2 = len - len1;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;
nFFT = 2 * len;

x_old = zeros(len1,1);
%n_old = zeros(len1,1);
xfinal = zeros(Nframes*len2, 1);

k = 1;
img = sqrt(-1);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

mu1 = 0.8;
mu2 = 0.6;
mu3 = 0.7;

threshold1 = 0.2;
%threshold1 = 0.3;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the clean speech
m1 = -18.4705;
m2 = -21.4016;
v1 = 18.0409;
v2 = 49.9072;
w1 = 0.7772;
w2 = 0.2228;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the high noise power
m3 = -46.2442;
m4 = -58.1222;
v3 = 32.9259;
v4 = 70.9205;
w3 = 0.7777;
w4 = 0.2223;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the GMM paramaters for the low noise power
m5 = -53.6041;
m6 = -45.0366;
v5 = 72.7698;
v6 = 27.0183;
w5 = 0.3750;
w6 = 0.6250;
w6 = 0.6250;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

count1 = 0;
count2 = 0;
count3 = 0;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    insign = f(n,:)';

    % Take Fourier transform of  frame
    spec = fft(insign, nFFT);

    % Compute the magnitude
    sig = abs(spec); 

    % Find the power spectrum
    %sig2 = spec .* conj(spec);
    sig2 = sig .^ 2;
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------
    
    % Assume the first frame is noise
    if (n == 1)
        noise_mu2 = sig2;
        final = noise_mu2;
    end
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------
    
%     % Compute the a posteriori SNR
%     gammak = min(sig2./noise_mu2, 40);  
% 
%     if n==1
%         ksi=aa+(1-aa)*max(gammak-1,0);
%     else
%         ksi=aa*Xk_prev./noise_mu2 + (1-aa)*max(gammak-1,0);     
% 
%         % decision-direct estimate of a priori SNR
%         ksi=max(ksi_min,ksi);  % limit ksi to -25 dB
%     end
% 
%     log_sigma_k = gammak.* ksi./ (1+ ksi)- log(1+ ksi); 
%     vad_decision = sum(log_sigma_k)/nFFT;    
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------

    % Noise only frame found
    %if (vad_decision < eta)
    if (pv(n) < 0.5)
        %mu = 0.9;
        mu = mu1;
        %noise_mu2 = mu * noise_mu2 + (1 - mu) * sig2;
        noise_mu2 = mu * sig2;
    else
        % Initialize the matrix
        noise_mu2_total = [];
        
        % Use 10 sub-frames 
        for subFrame = 1 : number_subFrames 
            % Define the input signal
            insign = f(n,:)';
            insign2 = insign((1 + ((subFrame-1)*(length(insign)/10))) : (subFrame * (length(insign)/10)));

            % Take Fourier transform of  frame
            spec_2 = fft(insign2, nFFT);

            % Compute the magnitude
            sig_2 = abs(spec_2); 

            % Find the power spectrum
            %sig2 = spec .* conj(spec);
            sig2 = sig_2 .^ 2;
            %-------------------------------------------------------------------------------------------------------------------------------
            %-------------------------------------------------------------------------------------------------------------------------------

            % Initialize the matrices
            power_p = [];
            power_p2 = [];

           % Hypothesis H_0: Transient Noise exists
           for freq_loop = 1 : length(sig2)
                z = 10*log10(sig2(freq_loop));

                %power_p = [power_p (9.4970 * (exp(- (z + 67.5990)^2 / (180.7520))))];
                pdfCompute = pdf_compute(z, m1, m2, m3, m4, v1, v2, v3, v4, w1, w2, w3, w4);
                power_p = [power_p pdfCompute];
            end

            % Hypothesis H_1: Background Noise exists
            for freq_loop = 1 : length(sig2)
                z = 10*log10(sig2(freq_loop));

                %power_p2 = [power_p2 (9.5066 * (exp(- (z + 66.8998)^2 / (180.3850))))];
                pdfCompute2 = pdf_compute(z, m1, m2, m5, m6, v1, v2, v5, v6, w1, w2, w5, w6);
                power_p2 = [power_p2 pdfCompute2];
            end
            %-------------------------------------------------------------------------------------------------------------------------------
            %-------------------------------------------------------------------------------------------------------------------------------

            % We use a discrepancy function
            % Define the variables 
            final_function = 0;
            observed_noisy_speech = sqrt(sig2);
            previous_noise_estimate = noise_mu;

            % Implement the discrepancy function
            final_function = sum(abs(observed_noisy_speech-previous_noise_estimate)) / sum(previous_noise_estimate);
            final_function = abs(final_function);
            discrepancy_measure = min(final_function, 1);

            if (discrepancy_measure > 0.8)
                % We increase "power_p" by 5%
                power_p = power_p * 1.05;

                % We decrease "power_p2" by 5%
                power_p2 = power_p2 * 0.95;
            else 
                % We increase "power_p2" by 5%
                power_p2 = power_p2 * 1.05;

                % We decrease "power_p" by 5%
                power_p = power_p * 0.95;
            end
            %-------------------------------------------------------------------------------------------------------------------------------
            %-------------------------------------------------------------------------------------------------------------------------------

            lamda_l_k_frameBin = power_p2 ./ power_p;

            % Use the geometric mean
            % find lamda_l_frame
            lamda_l_frame = 0;
            for u = 1 : length(sig2)
               lamda_l_frame = lamda_l_frame * lamda_l_k_frameBin(u); 
            end

            % find lamda_l_frame
            lamda_l_frame = lamda_l_frame^(1 / length(sig2));
            %-------------------------------------------------------------------------------------------------------------------------------
            %-------------------------------------------------------------------------------------------------------------------------------

            r_aPriori = 0.6 / 0.4; 

            % find final 
            final = ((r_aPriori * lamda_l_frame) / (1 + (r_aPriori * lamda_l_frame))) * (10*log10(final))' + (1 / (1 + (r_aPriori * lamda_l_frame))) * (10*log10(sig2))'; 
            final = final';
            final = power(10, (final / 10));        
            %-------------------------------------------------------------------------------------------------------------------------------
            %-------------------------------------------------------------------------------------------------------------------------------

            % We use the discrepancy function
            if (discrepancy_measure > 0.8)
                % Update the noise power spectrum
                %mu = 0.2;
                mu = mu2;
                
                %noise_mu2 = mu * noise_mu2 + (1 - mu) * final;
                noise_mu2 = (1-mu) * final;
            else
                % Update the noise power spectrum
                %mu = 0.9;
                mu = mu3;
                
                %noise_mu2 = mu * noise_mu2 + (1 - mu) * final;
                noise_mu2 = (1-mu) * final;
            end

            % Define the noise power of the sub-frames
            noise_mu2_total = [noise_mu2_total; noise_mu2'];
        end

        % Define the final noise power estimate
        % Use the average of the sub-frames
        noise_mu2 = mean(noise_mu2_total, 1);
        noise_mu2 = noise_mu2';
    end
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------

    % Define "noise_mu"
    noise_mu = sqrt(noise_mu2);

    % Store the noise power
    total_noise_power_frame = [total_noise_power_frame; noise_mu2'];
    %-------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------

%     vk=ksi.*gammak./(1+ksi);
%     j0 = besseli(0,vk/2);
%     j1 = besseli(1,vk/2);
% 
%     C=exp(-0.5*vk);
%     A=((c*(vk.^0.5)).*C)./gammak;
%     B=(1+vk).*j0+vk.*j1;
%     hw=A.*B;
% 
%     sig=sig.*hw;
% 
%     % save for estimation of a priori SNR in next frame
%     Xk_prev=sig.^2;  
%     %-------------------------------------------------------------------------------------------------------------------------------
%     %-------------------------------------------------------------------------------------------------------------------------------
% 
%     xi_w= ifft( sig .* exp(img*angle(spec)),nFFT);
% 
%     xi_w= real( xi_w);
%     %-------------------------------------------------------------------------------------------------------------------------------
%     %-------------------------------------------------------------------------------------------------------------------------------
% 
%     % --- Overlap and add ---------------
%     xfinal(k:k+ len2-1)= x_old+ xi_w(1:len1);
%     x_old= xi_w(len1+ 1: len);
% 
%     % --- Overlap and add ---------------
%     ni_time = real(ifft(noise_mu));
%     ni_time_final(k:k+len2-1) = n_old + ni_time(1:len1);
%     n_old = ni_time(1+len1:len);
% 
%     % Update k
%     k = k + len2;
end
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% We use MSE
MSE = 0;

% Check the result of the algorittm
% We use: "total_noise_power_frame"
for i = 1 : size(total_noise_power_frame, 1)
    for k = 1: size(total_noise_power_frame, 2)
        % Compute the MSE
        MSE = MSE + ((total_noise_power_frame(i,k) - true_total_noise_power_frame(i,k))^2 / true_total_noise_power_frame(i,k));
    end
end

% Find the final value of the MSE
MSE = MSE / (size(total_noise_power_frame, 1)*size(total_noise_power_frame, 2));
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
















%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

filename2 = 'SA2.WAV';
filename = 'buccaneer2';
snr = 0; 
integrateToSS = 0; 
integrateToMMSE = 1;
useGMM_logPower = 0;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal, returnFromMMSE] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, 0, 0);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------























%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 10; 
integrateToSS = 1; 
integrateToMMSE = 0;
plotGMMs = 0;
plotPowerTime = 1;
useGMM_logPower = 2;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
hold on;


% Read the input file
filename2 = 'SA1.WAV'; 
fs = 16000;
s = readsph(fullfile(filename2));

% Read the noise file
filename = 'destroyerengine'; 
[s2, fs2] = readwav(fullfile(filename));

snr = 10;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'dopEk',s2,fs2); 

% We use this for noise only: "z2(:,2)"
z2 = v_addnoise(s,fs,snr,'dxopEk',s2,fs2);
% We use "soundsc(z2(:,2),fs)" for noise only
% Read the input file
filename2 = 'SA1.WAV'; fs = 16000;
s = readsph(fullfile(filename2));

% Read the noise file
filename = 'destroyerengine'; 
[s2, fs2] = readwav(fullfile(filename));

% Define the SNR
snr = 0;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = [z; v_addnoise(s,fs,snr,'dopEk',s2,fs2)]; 

% We use this for noise only: "z2(:,2)"
z2 = [z2; v_addnoise(s,fs,snr,'dxopEk',s2,fs2)];

method = 'mcra';
specsub_ns22(z, method);
legend('Proposed','True','MCRA');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Call the main function 
filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 15; 
integrateToSS = 1; 
integrateToMMSE = 0;
plotGMMs = 0;
plotPowerTime = 0;
useGMM_logPower = 0;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Check the result from SS
%soundsc(result_fromSS,16000);

% We can listen to the noise signal. But, this does not really make much sense.
%soundsc(result_fromSS_noise,16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the clean speech signal
fs = 16000;
writewav(clean_speech_signal, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final1 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Call the main function 
filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 10; 
integrateToSS = 1; 
integrateToMMSE = 0;
plotGMMs = 0;
plotPowerTime = 0;
useGMM_logPower = 0;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Check the result from SS
%soundsc(result_fromSS,16000);

% We can listen to the noise signal. But, this does not really make much sense.
%soundsc(result_fromSS_noise,16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the clean speech signal
fs = 16000;
writewav(clean_speech_signal, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final2 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Call the main function 
filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 5; 
integrateToSS = 1; 
integrateToMMSE = 0;
plotGMMs = 0;
plotPowerTime = 0;
useGMM_logPower = 0;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Check the result from SS
%soundsc(result_fromSS,16000);

% We can listen to the noise signal. But, this does not really make much sense.
%soundsc(result_fromSS_noise,16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the clean speech signal
fs = 16000;
writewav(clean_speech_signal, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final3 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Call the main function 
filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 0; 
integrateToSS = 1; 
integrateToMMSE = 0;
plotGMMs = 0;
plotPowerTime = 0;
useGMM_logPower = 0;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Check the result from SS
%soundsc(result_fromSS,16000);

% We can listen to the noise signal. But, this does not really make much sense.
%soundsc(result_fromSS_noise,16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the clean speech signal
fs = 16000;
writewav(clean_speech_signal, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final4 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
T = table(final1, final2, final3, final4);
T.Properties.VariableNames = {'TestResult1', 'TestResult2', 'TestResult3', 'TestResult4'};

% print the results
disp(['The following table shows the results that we just calculated.']);
disp(T);
disp(['   ']); disp(['   ']);

% set the figure data
f = figure; set(gcf,'Color','w');
set(f, 'Position', [500 400 300 100]);
dat =  {T{1,1}, T{1,2}, T{1,3}, T{1,4};...
        T{2,1}, T{2,2}, T{2,3}, T{2,4};...   
        T{3,1}, T{3,2}, T{3,3}, T{3,4};...
        T{4,1}, T{4,2}, T{4,3}, T{4,4};...
        T{5,1}, T{5,2}, T{5,3}, T{5,4};...
        T{6,1}, T{6,2}, T{6,3}, T{6,4};...
        T{7,1}, T{7,2}, T{7,3}, T{7,4};...
        T{8,1}, T{8,2}, T{8,3}, T{8,4};...
        T{9,1}, T{9,2}, T{9,3}, T{9,4};...
        T{10,1}, T{10,2}, T{10,3}, T{10,4};...
        T{11,1}, T{11,2}, T{11,3}, T{11,4};...
        T{12,1}, T{12,2}, T{12,3}, T{12,4};...
        T{13,1}, T{13,2}, T{13,3}, T{13,4};...
        T{14,1}, T{14,2}, T{14,3}, T{14,4};...
        T{15,1}, T{15,2}, T{15,3}, T{15,4};...
        T{16,1}, T{16,2}, T{16,3}, T{16,4};...
        T{17,1}, T{17,2}, T{17,3}, T{17,4}; ...
        T{18,1}, T{18,2}, T{18,3}, T{18,4}; ...
        T{19,1}, T{19,2}, T{19,3}, T{19,4}; ...
        T{20,1}, T{20,2}, T{20,3}, T{20,4}; ...
        T{21,1}, T{21,2}, T{21,3}, T{21,4};};

% set the table as a figure
columnname =   {'Results for 15 dB SNR', 'Results for 10 dB SNR', 'Results for 5 dB SNR', 'Results for 0 dB SNR'};
rnames = {'SNR mean', 'Segmental SNR', 'WSS mean', 'LLR mean', 'IS mean', 'CEP mean', ...
    'fw SNR segmental', 'SIG' 'BAK', 'OVL', 'SIG2' 'BAK2', 'OVL2', 'raw PESQ', 'Csig', 'Cbak', 'Covl', 'NCM' , 'CSh', 'CSm', 'CSl'};
columnformat = {'numeric'}; 
t = uitable('Units','normalized','Position',...
            [0.05 0.05 0.755 0.87], 'Data', dat,... 
            'ColumnName', columnname,...
            'RowName', rnames, ...
            'ColumnFormat', columnformat); 
        
%set the title
%set(gca,'FontSize',15); 
title('Table: Proposed Algorithm with SS. Only one speaker and one noise is tested.');        
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------












%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Call the main function 
filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 5; 
integrateToSS = 1; 
integrateToMMSE = 0;
plotGMMs = 0;
plotPowerTime = 1;
useGMM_logPower = 2;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal, a] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
result_fromSS = a;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Check the result from SS
%soundsc(result_fromSS,16000);

% We can listen to the noise signal. But, this does not really make much sense.
%soundsc(result_fromSS_noise,16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the clean speech signal
fs = 16000;
writewav(clean_speech_signal, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final1 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Call the main function 
filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 10; 
integrateToSS = 0; 
integrateToMMSE = 1;
plotGMMs = 0;
plotPowerTime = 0;
useGMM_logPower = 0;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal, a] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
result_fromSS = a;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Check the result from SS
%soundsc(result_fromSS,16000);

% We can listen to the noise signal. But, this does not really make much sense.
%soundsc(result_fromSS_noise,16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the clean speech signal
fs = 16000;
writewav(clean_speech_signal, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final2 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Call the main function 
filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 5; 
integrateToSS = 0; 
integrateToMMSE = 1;
plotGMMs = 0;
plotPowerTime = 0;
useGMM_logPower = 0;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal, a] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
result_fromSS = a;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Check the result from SS
%soundsc(result_fromSS,16000);

% We can listen to the noise signal. But, this does not really make much sense.
%soundsc(result_fromSS_noise,16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the clean speech signal
fs = 16000;
writewav(clean_speech_signal, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final3 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Call the main function 
filename2 = 'SA1.WAV';
filename = 'destroyerengine';
snr = 0; 
integrateToSS = 0; 
integrateToMMSE = 1;
plotGMMs = 0;
plotPowerTime = 0;
useGMM_logPower = 0;
[total_noise_power_frame, MSE, time_T, result_fromSS, result_fromSS_noise, clean_speech_signal, a] = proposed_noiseEstimation(filename2, filename, snr, useGMM_logPower, integrateToSS, integrateToMMSE, plotGMMs, plotPowerTime);
result_fromSS = a;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Check the result from SS
%soundsc(result_fromSS,16000);

% We can listen to the noise signal. But, this does not really make much sense.
%soundsc(result_fromSS_noise,16000);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Store the clean speech signal
fs = 16000;
writewav(clean_speech_signal, fs, 'clean_speech_signal_main');

% Store the enhanced speech signal
writewav(result_fromSS, fs, 'result_fromSS_speech_signal_main');

% Add the datapath for Quality Measures
addpath ./Functions_for_Noise_Estimation\quality

% Use comp_snr
% Find the Overall and Segmental SNR
[snr_mean, segsnr_mean] = comp_snr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_wss
% Find the Weighted-Spectral Slope 
wss_mean = comp_wss('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_llr
% Find the Likelihood Ratio measure
llr_mean = comp_llr('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_is
% Find the Itakura-Saito measure
is_mean = comp_is('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_cep
% Find the cepstral distance measure
cep_mean = comp_cep('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use comp_fwseg
% Find the Frequency Weighted Semgmental SNR
fwSNRseg = comp_fwseg('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use comp_fwseg_variant
% Find the Frequency-variant fwSNRseg measure
[SIG, BAK, OVL] = comp_fwseg_variant('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'SIG' is the predicted rating of speech distortion,
%           'BAK' is the predicted rating of background noise distortion,
%           'OVL' is the predicted rating of overall quality.

% Use comp_fwseg_mars
% Find the Frequency variant fwSNRseg measure based on MARS analysis
[SIG2, BAK2, OVL2] = comp_fwseg_mars('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Use the PESQ measure
pesq_val = pesq('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use composite
% Find a composite emasure (using certain already found metrics)
[Csig, Cbak, Covl] = composite('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');
%	where   'Csig' is the predicted rating of speech distortion,
%           'Cbak' is the predicted rating of background noise distortion,
%           'Covl' is the predicted rating of overall quality.
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Add the datapath for Intelligibility Measures
addpath ./Functions_for_Noise_Estimation\intelligibility

% Use the NCM measure
ncm = NCM('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the CSII measure
[csh, csm, csl] = CSII('clean_speech_signal_main.wav', 'result_fromSS_speech_signal_main.wav');

% Use the SII measure
% EXAMPLE
% sp=[40 45 50 24 56 60 55 55 52 48 50 51 55 67 76 67 56 31];
% ns=[30 50 60 20 60 50 70 45 80 40 60 20 60 22 55 50 67 40];
% M= 5;
% sv = SII (sp,ns, M);
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Define the raw PESQ value
raw_PESQ = pesq_val;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
final = [snr_mean, segsnr_mean, wss_mean, llr_mean, is_mean, cep_mean, ...
    fwSNRseg, SIG, BAK, OVL, SIG2, BAK2, OVL2, raw_PESQ, Csig, Cbak, Covl, ...
    ncm, csh, csm, csl];
final4 = final';
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% create a table with the values calculated
T = table(final1, final2, final3, final4);
T.Properties.VariableNames = {'TestResult1', 'TestResult2', 'TestResult3', 'TestResult4'};

% print the results
disp(['The following table shows the results that we just calculated.']);
disp(T);
disp(['   ']); disp(['   ']);

% set the figure data
f = figure; set(gcf,'Color','w');
set(f, 'Position', [500 400 300 100]);
dat =  {T{1,1}, T{1,2}, T{1,3}, T{1,4};...
        T{2,1}, T{2,2}, T{2,3}, T{2,4};...   
        T{3,1}, T{3,2}, T{3,3}, T{3,4};...
        T{4,1}, T{4,2}, T{4,3}, T{4,4};...
        T{5,1}, T{5,2}, T{5,3}, T{5,4};...
        T{6,1}, T{6,2}, T{6,3}, T{6,4};...
        T{7,1}, T{7,2}, T{7,3}, T{7,4};...
        T{8,1}, T{8,2}, T{8,3}, T{8,4};...
        T{9,1}, T{9,2}, T{9,3}, T{9,4};...
        T{10,1}, T{10,2}, T{10,3}, T{10,4};...
        T{11,1}, T{11,2}, T{11,3}, T{11,4};...
        T{12,1}, T{12,2}, T{12,3}, T{12,4};...
        T{13,1}, T{13,2}, T{13,3}, T{13,4};...
        T{14,1}, T{14,2}, T{14,3}, T{14,4};...
        T{15,1}, T{15,2}, T{15,3}, T{15,4};...
        T{16,1}, T{16,2}, T{16,3}, T{16,4};...
        T{17,1}, T{17,2}, T{17,3}, T{17,4}; ...
        T{18,1}, T{18,2}, T{18,3}, T{18,4}; ...
        T{19,1}, T{19,2}, T{19,3}, T{19,4}; ...
        T{20,1}, T{20,2}, T{20,3}, T{20,4}; ...
        T{21,1}, T{21,2}, T{21,3}, T{21,4};};

% set the table as a figure
columnname =   {'Results for 15 dB SNR', 'Results for 10 dB SNR', 'Results for 5 dB SNR', 'Results for 0 dB SNR'};
rnames = {'SNR mean', 'Segmental SNR', 'WSS mean', 'LLR mean', 'IS mean', 'CEP mean', ...
    'fw SNR segmental', 'SIG' 'BAK', 'OVL', 'SIG2' 'BAK2', 'OVL2', 'raw PESQ', 'Csig', 'Cbak', 'Covl', 'NCM' , 'CSh', 'CSm', 'CSl'};
columnformat = {'numeric'}; 
t = uitable('Units','normalized','Position',...
            [0.05 0.05 0.755 0.87], 'Data', dat,... 
            'ColumnName', columnname,...
            'RowName', rnames, ...
            'ColumnFormat', columnformat); 
        
%set the title
%set(gca,'FontSize',15); 
title('Table: Proposed Algorithm with MMSE. Only one speaker and one noise is tested.');        
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% tracking delay

filename = 'sp01_babble_sn15';
[s,fs] = readwav(fullfile(filename));

z = s;
%z = [z; z];

start_point1 = length(s);

filename = 'sp01_16kHz';
[s2,fs2] = readwav(fullfile(filename));

rrrqae = s - s2;
%rrrqae = [rrrqae; rrrqae];

filename = 'sp01_babble_sn0';
[s,fs] = readwav(fullfile(filename));

z = [z; s; s];

filename = 'sp01_16kHz';
[s2,fs2] = readwav(fullfile(filename));

rrrqae2 = s - s2;
rrrqae = [rrrqae; rrrqae2; rrrqae2];

method = 'mcra2';
[mcra2_noise_ps] = new_new2_noise_parameters(z, method);

true_ppx = true25_noise_parameters_my_proposed_algorithm( rrrqae );
true_ppx = 10*log10(true_ppx);


value_main = [];
for tr = 1 : length(mcra2_noise_ps(1,:))
    new_value = mean(mcra2_noise_ps(end-50:end,tr));

    for kk = length(mcra2_noise_ps(:,1))/2 : length(mcra2_noise_ps(:,1))
        if (mcra2_noise_ps(kk,tr) > 0.9*new_value)
            value_main = [value_main kk];
        end
    end
end

value_main = round(mean(value_main));

method = 'conn_freq';
[conn_freq_noise_ps] = new_new2_noise_parameters(z, method);

method = 'imcra';
[imcra_noise_ps] = new_new2_noise_parameters(z, method);

method = 'mcra';
[mcra_noise_ps] = new_new2_noise_parameters(z, method);

method = 'martin';
[martin_noise_ps] = new_new2_noise_parameters(z, method);

method = 'doblinger';
[doblinger_noise_ps] = new_new2_noise_parameters(z, method);

method = 'hirsch';
[hirsch_noise_ps] = new_new2_noise_parameters(z, method);

true_ppx = true25_noise_parameters_my_proposed_algorithm( rrrqae );

[yp ~] = size(hirsch_noise_ps);

atasgas = [];
atasgas2 = [];
atasgas3 = [];
atasgas4 = [];
atasgas5 = [];
atasgas6 = [];
atasgas7 = [];
atasgas8 = [];

for i = 1 : yp
    if (i == 1) 
        var = 0.030;
    end
    
    atasgas = [atasgas 10*log10(mcra2_noise_ps(i,20))]; 
    atasgas2 = [atasgas2 10*log10(conn_freq_noise_ps(i,20))]; 
    atasgas3 = [atasgas3 10*log10(imcra_noise_ps(i,20))]; 
    atasgas4 = [atasgas4 10*log10(mcra_noise_ps(i,20))]; 
    atasgas5 = [atasgas5 10*log10(martin_noise_ps(i,20))]; 
    atasgas6 = [atasgas6 10*log10(doblinger_noise_ps(i,20))]; 
    atasgas7 = [atasgas7 10*log10(hirsch_noise_ps(i,20))]; 
    atasgas8 = [atasgas8 10*log10(true_ppx(i,20))]; 

    var = var + (0.030/2);
end

figure; set(gcf,'Color','w');

set(gca,'FontSize',15);
plot(0.020:(0.030/2):var-(0.030/2), atasgas);
hold on;

plot(0.020:(0.030/2):var-(0.030/2), atasgas2, 'r');
hold on;

% plot(0.020:(0.030/2):var-(0.030/2), atasgas3, 'g');
% hold on;

plot(0.020:(0.030/2):var-(0.030/2), atasgas4, 'g');
hold on;

plot(0.020:(0.030/2):var-(0.030/2), atasgas5, 'm');
hold on;

plot(0.020:(0.030/2):var-(0.030/2), atasgas6, 'c');
hold on;

% plot(0.020:(0.030/2):var-(0.030/2), atasgas7, 'k');
% hold on;

plot(0.020:(0.030/2):var-(0.030/2), atasgas8, 'k');
hold off;

set(gca,'FontSize',15);
title('Noise Estimation Algorithms. 30 ms with 15 ms overlap. Power (in dB) vs Time (s) for one frequency bin.');

xlabel('Time in seconds (s)');
ylabel('Power in dB');

legend(' MCRA2',' Connected TF', ' MCRA', ' Martin MS', ' Doblinger', ' True', 'Location', 'SouthEast');

grid on; grid minor;
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

































%%
% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% Read the noise file
[s2, fs2, bits] = readwav('babble.wav');

% Define the SNR
snr = 15;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

% s is the clean speech signal
% z is the noisy speech signal

%time = (1/fs) : (1/fs) : ((length(s)*(1/fs))-(1/fs));
time = (1/fs) : (1/fs) : ((length(s)*(1/fs)));

% plot(time, s);

figure; 
set(gcf,'Color','w');
subplot(2,1,1);
set(gca,'FontSize',15);

plot(time, s);

set(gca,'FontSize',15);
title('Clean speech signal.');
xlabel('Time');
ylabel('Amplitude');

set(gcf,'Color','w');
subplot(2,1,2);

plot(time, z);

set(gca,'FontSize',15);
title('Noisy speech signal.');
xlabel('Time');
ylabel('Amplitude');

% print the time 
time(end)

% we have 3 seconds 
% we need 10 seconds 
z = [z; z; z; z];
s = [s; s; s; s];

time = (1/fs) : (1/fs) : ((length(s)*(1/fs)));

figure; 
set(gcf,'Color','w');
subplot(2,1,1);
set(gca,'FontSize',15);

plot(time, s);

set(gca,'FontSize',15);
title('Clean speech signal.');
xlabel('Time');
ylabel('Amplitude');

set(gcf,'Color','w');
subplot(2,1,2);

plot(time, z);

set(gca,'FontSize',15);
title('Noisy speech signal.');
xlabel('Time');
ylabel('Amplitude');

% print the time 
time(end)

x = find(time>10,1);

%x = x + 1;

time(x)

time = time(1:x);

%time = [0 time];

s = s(1:x);
z = z(1:x);

figure; 
set(gcf,'Color','w');
subplot(2,1,1);
set(gca,'FontSize',15);

plot(time, s);

set(gca,'FontSize',15);
title('Clean speech signal.');
xlabel('Time');
ylabel('Amplitude');

set(gcf,'Color','w');
subplot(2,1,2);

plot(time, z);

set(gca,'FontSize',15);
title('Noisy speech signal.');
xlabel('Time');
ylabel('Amplitude');

% print the time 
time(end)

% we have thus the input signal z
% z is 10 seconds
% z is the noisy speech signal

total_observations = [];

% Define the frame size in samples
% len = txinc_samples*10;
len = 0.090 * fs;

% window overlap in percent of frame size
PERC = (80/90) * 100;

len1 = floor(len * PERC/100);
len2 = len - len1;

x = z;

% Define the parameters
Nframes = floor(length(x)/len2)-1;

nFFT = 2 * len;

k = 1;
img = sqrt(-1);

% Define window
win = hamming(len);

% Nframes = 92;

% Main for loop
for n = 1 : Nframes
    
    if (k+len-1 <= length(x))
        % Define the input signal
        signal_part = x(k:k+len-1);
        insign = signal_part .* win;

        % Take Fourier transform of  frame
        spec = fft(insign, nFFT);

        sig = abs(spec);
        sig2 = sig.^2;

        phase_t = phase(spec);

        total_observations = [total_observations; phase_t'];

        k = k+len2;
    end
end








asfasdfasdfas







% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


%%
% KF for pv

% Things to do: 
% 1) EM algorithm VS LS/MMSE for training the KF. We use only LS/MMSE here. We need to use the EM algorithm.

% 2) Use future frames. Use KF for smoothing/hindsight and not for filtering. In general, filtering VS prediction VS smoothing.

% 3) Now, pv is the frame voiced probability. We need to try pv(k,l), we need to try fr.bin and frame voiced probability.

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we store the true values
true_pv = pv;

% we use noise 

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% Read the noise file
[s2, fs2, bits] = readwav('babble.wav');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we find error
error_pv = pv - true_pv;

% fit zero mean gaussian
[m,v,w] = gaussmix(error_pv, [], [], 1);
[m2,v2,w2] = gaussmix(error_pv-m, [], [], 1);

% we use v2
v2

% clear the other temp variables
clear m;
clear v;
clear w;
clear m2;
clear w2;

% this was for one sentence
% we need to repeat with many sentences

% we initialize the counter
counter = 0;

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    total_error_pv = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 5;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        % we find error
        error_pv = pv - true_pv;
        
        total_error_pv = [total_error_pv; error_pv'];
        
        total_total_error_pv{counter} = total_error_pv;
        
    end
end

% this was for 5 dB SNR 
% we now use 0 dB SNR 

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    total_error_pv = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 0;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        % we find error
        error_pv = pv - true_pv;
        
        total_error_pv = [total_error_pv; error_pv'];
        
        total_total_error_pv{counter} = total_error_pv;
        
    end
end

% this was for 0 dB SNR 
% we now use 10 dB SNR 

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    total_error_pv = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 10;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        % we find error
        error_pv = pv - true_pv;
        
        total_error_pv = [total_error_pv; error_pv'];
        
        total_total_error_pv{counter} = total_error_pv;
        
    end
end

% this was for 0 dB SNR 

% we now need one array with all the values

% we initialize the new array
errors_in_pv = [];

for i = 1 : counter
    temp_array = total_total_error_pv{i};
    
    for k = 1 : size(temp_array,1)
        error_in_pv = [errors_in_pv temp_array(k,:)];
    
    end
end

% we use: error_in_pv

error_in_pv = error_in_pv';

% fit zero mean gaussian
[m,v,w] = gaussmix(error_in_pv, [], [], 1);
[m2,v2,w2] = gaussmix(error_in_pv-m, [], [], 1);

% we use v2
v2

% v2 is 0.0942
% we use this value so as to not wait for training
%v2 = 0.0942;

% clear the other temp variables
clear m;
clear v;
clear w;
clear m2;
clear w2;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now fit the transition model 

tr_model = [];

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    tr_model = [tr_model; diff(pv)];

end

% fit zero mean gaussian
[m,v,w] = gaussmix(tr_model, [], [], 1);
[m2,v2_x,w2] = gaussmix(tr_model-m, [], [], 1);

% we use v2_x
v2_x

% v2_x is 0.0198
% we use this value so as to not wait for training
%v2_x = 0.0198;

% clear the other temp variables
clear m;
clear v;
clear w;
clear m2;
clear w2;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now use Kalman filter (KF)

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

reference_true_pv = pv;

% Read the noise file
[s2, fs2] = readwav('babble.wav');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

predicted_x_m = [];
predicted_P_m = [];

for loop = 1 : length(pv)
    % we initialize x_m and P_m
    if (loop == 1)
        x_m = 0;
        P_m = 10.2;
        
    end
    
    x_m = (((P_m+v2_x)*pv(loop)) + (v2*x_m)) / (P_m+v2_x+v2);
    
    P_m = ((P_m+v2_x) * v2) / (P_m+v2_x+v2);
    
    % store the results
    predicted_x_m = [predicted_x_m x_m];
    predicted_P_m = [predicted_P_m P_m];

end

predicted_x_m = predicted_x_m';
predicted_P_m = predicted_P_m';

figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);

plot(predicted_x_m);
hold on;
plot(pv, 'g');
hold on;
plot(reference_true_pv, 'r');
hold off;

set(gca,'FontSize',10);
title('pv. Probability of frame being voiced.');
xlabel('Time in frames');
ylabel('Probability pv');

legend(' new pv', ' PEFAC pv', ' true pv');

% ------------------------------------------------------------------------

% treat the problem as a classification problem
% we use 0 and 1 only

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) > 0.5)
       reference_true_pv(i) = 1;
   else
       reference_true_pv(i) = 0;
   end
end

for i = 1 : length(pv)
   if (pv(i) > 0.5)
       pv(i) = 1;
   else
       pv(i) = 0;
   end
end

for i = 1 : length(predicted_x_m)
   if (predicted_x_m(i) > 0.5)
       predicted_x_m(i) = 1;
   else
       predicted_x_m(i) = 0;
   end
end

% we find accuracy for pv
pv_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == pv(i))
       pv_accuracy = pv_accuracy + 1;
   end
end

pv_accuracy = pv_accuracy / length(reference_true_pv);

% we find accuracy for x_m
x_m_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == predicted_x_m(i))
       x_m_accuracy = x_m_accuracy + 1;
   end
end

x_m_accuracy = x_m_accuracy / length(reference_true_pv);

% we now use pv
% we compare pv and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

pv_true_positives = 0;
pv_true_negatives = 0;

pv_false_positives = 0;
pv_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && pv(i) == 1)
        pv_true_positives = pv_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 0)
        pv_true_negatives = pv_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && pv(i) == 0)
        pv_false_positives = pv_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 1)
        pv_false_negatives = pv_false_negatives + 1;
    end
end
    
% find precision
pv_precision = (pv_true_positives) / (pv_true_positives+pv_false_positives);

% find recall
pv_recall = (pv_true_positives) / (pv_true_positives+pv_false_negatives);

% we now do the same for x_m

% we now use x_m
% we compare x_m and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

x_m_true_positives = 0;
x_m_true_negatives = 0;

x_m_false_positives = 0;
x_m_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && predicted_x_m(i) == 1)
        x_m_true_positives = x_m_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 0)
        x_m_true_negatives = x_m_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && predicted_x_m(i) == 0)
        x_m_false_positives = x_m_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 1)
        x_m_false_negatives = x_m_false_negatives + 1;
    end
end

% find precision
x_m_precision = (x_m_true_positives) / (x_m_true_positives+x_m_false_positives);

% find recall
x_m_recall = (x_m_true_positives) / (x_m_true_positives+x_m_false_negatives);

% we thus have 3 metrics for the two methods

% we have: pv_accuracy, pv_precision, pv_recall 
% we have: x_m_accuracy, x_m_precision, x_m_recall

% print the metrics
disp(['The pv PEFAC has accuracy: ', num2str(pv_accuracy)]);
disp(['The pv PEFAC has precision: ', num2str(pv_precision)]);
disp(['The pv PEFAC has recall: ', num2str(pv_recall)]);

disp([' ']);
disp([' ']);

% print the metrics
disp(['The new pv method has accuracy: ', num2str(x_m_accuracy)]);
disp(['The new pv method has precision: ', num2str(x_m_precision)]);
disp(['The new pv method has recall: ', num2str(x_m_recall)]);

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now repeat the KF

% we use a different training method
% we use LS or MMSE training for the KF

% we use non-iterative approach 

% we need to define total_total_observations 
% we need to define the states: total_states

% time increases in the column direction 

counter = 0;

states1 = [];
observations2 = [];

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    states2 = [];
    observations2 = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 5;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        states2 = [states2; true_pv'];
        observations2 = [observations2; pv'];
        
        states1{counter} = states2;
        observations1{counter} = observations2;
        
    end
end

% this was for SNR 5 dB
% we now use SNR 0 dB

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    states2 = [];
    observations2 = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 0;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        states2 = [states2; true_pv'];
        observations2 = [observations2; pv'];
        
        states1{counter} = states2;
        observations1{counter} = observations2;
        
    end
end

% this was for SNR 0 dB
% we now use SNR 10 dB

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    states2 = [];
    observations2 = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 10;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        states2 = [states2; true_pv'];
        observations2 = [observations2; pv'];
        
        states1{counter} = states2;
        observations1{counter} = observations2;
        
    end
end

% SNR 5, 0 and 10 were used
% the training with SNR 5,0 and 10 has finished

% we need to define total_total_observations 
% we need to define the states: total_states

total_states = [];
total_total_observations = [];

temp_array = [];

for i = 1 : counter
    temp_array = states1{i};
    temp_array2 = observations1{i};
    
    for k = 1 : size(temp_array,1)
        for kk = 1 : size(temp_array,2)
            total_states = [total_states temp_array(k,kk)];
            total_total_observations = [total_total_observations temp_array2(k,kk)];
        
        end
    end
end

% we defined total_total_observations 
% we defined the states: total_states

% we now need to train the KF

% ------------------------------------------------------------------------

% KF -- train parameters -- assume constant

% we use LS training (or MMSE training and NOT the EM algorithm)

% In KF, there is no access to the correct data when the algorithm is running in real time.
% In KF, the error is between the measurement and the Kalman prediction.
% It can be stated that the KF resembles a low pass filter.
% but its transfer characteristic is nonlinear, and the cut-off frequency shifts.

% In KF, there are 2 main things to do in each time step. 
% The first thing is to calculate the apriori state estimate, and the measurement update.
% The second step is to calculate the gain, corresponding covariance, and then the posterior state.

% we use M time steps
M = size(total_states,2);

for loop = 1 : M
    temp = total_states(:,loop);
    
    X(:,loop) = temp';
end

X1 = X(:, 1:(end-1));
X2 = X(:, 2:end);

clear temp;
C = 1;

for loop = 1 : M
    temp = zeros(1,C);

    % temp = observations;
    temp = total_total_observations(:,loop);

    Z(:,loop) = temp';
end

%A = X2 * X1' * (X1 * X1')^(-1);
A = X2 * X1' * pinv(X1 * X1');

%H = Z * X' * (X * X')^(-1);
H = Z * X' * pinv(X * X');

W = (X2 - (A * X1)) * (X2 - (A * X1))' / (M-1);

Q = (Z - (H * X)) * (Z - (H * X))' / M; 

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% do the KF

% A is the same as the F matrix
% W is the same as the Q matrix

F = A;
Q1 = Q;
Q = W;

% H is the same as H
% Q is the same as the R matrix

R = Q1;

% we use noisy speech signal
% we implement the KF filter

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

reference_true_pv = pv;

% Read the noise file
[s2, fs2] = readwav('babble.wav');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

predicted_x_m = [];
predicted_P_m = [];

for loop = 1 : length(pv)
    % we initialize x_m and delta_delta_phi and P_p
    if (loop == 1)
        x_m = 0.5;
        P_m = 1000.2;
    end

    % define l 
    % l is the observation
    l = pv(loop);
    
    K = ((F*P_m*F') + Q) * H' * ((H * ((F*P_m*F') + Q) * H') + R)^(-1);
    
    x_p = (F*x_m) + (K * (l - (H * F * x_m)));
    
    P_p = (1-K) * ((F * P_m * F') + Q);
    
    % constrain x_p
    if (x_p < 0)
       x_p = 0.01; 
    end
    
    % constrain x_p
    if (x_p > 1)
       x_p = 0.99; 
    end
    
    % store the results
    predicted_x_m = [predicted_x_m x_p];
    predicted_P_m = [predicted_P_m P_p];

    % prepare the next iteration
    x_m = x_p;
    P_m = P_p;

end

predicted_x_m = predicted_x_m';
predicted_P_m = predicted_P_m';

figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);

plot(predicted_x_m);
hold on;
plot(pv, 'g');
hold on;
plot(reference_true_pv, 'r');
hold off;

set(gca,'FontSize',10);
title('pv. Probability of frame being voiced.');
xlabel('Time in frames');
ylabel('Probability pv');

legend(' new pv', ' PEFAC pv', ' true pv');

% ------------------------------------------------------------------------

% treat the problem as a classification problem
% we use 0 and 1 only

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) > 0.5)
       reference_true_pv(i) = 1;
   else
       reference_true_pv(i) = 0;
   end
end

for i = 1 : length(pv)
   if (pv(i) > 0.5)
       pv(i) = 1;
   else
       pv(i) = 0;
   end
end

for i = 1 : length(predicted_x_m)
   if (predicted_x_m(i) > 0.5)
       predicted_x_m(i) = 1;
   else
       predicted_x_m(i) = 0;
   end
end

% we find accuracy for pv
pv_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == pv(i))
       pv_accuracy = pv_accuracy + 1;
   end
end

pv_accuracy = pv_accuracy / length(reference_true_pv);

% we find accuracy for x_m
x_m_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == predicted_x_m(i))
       x_m_accuracy = x_m_accuracy + 1;
   end
end

x_m_accuracy = x_m_accuracy / length(reference_true_pv);

% we now use pv
% we compare pv and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

pv_true_positives = 0;
pv_true_negatives = 0;

pv_false_positives = 0;
pv_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && pv(i) == 1)
        pv_true_positives = pv_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 0)
        pv_true_negatives = pv_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && pv(i) == 0)
        pv_false_positives = pv_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 1)
        pv_false_negatives = pv_false_negatives + 1;
    end
end
    
% find precision
pv_precision = (pv_true_positives) / (pv_true_positives+pv_false_positives);

% find recall
pv_recall = (pv_true_positives) / (pv_true_positives+pv_false_negatives);

% we now do the same for x_m

% we now use x_m
% we compare x_m and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

x_m_true_positives = 0;
x_m_true_negatives = 0;

x_m_false_positives = 0;
x_m_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && predicted_x_m(i) == 1)
        x_m_true_positives = x_m_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 0)
        x_m_true_negatives = x_m_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && predicted_x_m(i) == 0)
        x_m_false_positives = x_m_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 1)
        x_m_false_negatives = x_m_false_negatives + 1;
    end
end

% find precision
x_m_precision = (x_m_true_positives) / (x_m_true_positives+x_m_false_positives);

% find recall
x_m_recall = (x_m_true_positives) / (x_m_true_positives+x_m_false_negatives);

% we thus have 3 metrics for the two methods

% we have: pv_accuracy, pv_precision, pv_recall 
% we have: x_m_accuracy, x_m_precision, x_m_recall

disp([' ']);
disp([' ']);

disp([' ']);
disp([' ']);

% print the metrics
disp(['The pv PEFAC has accuracy: ', num2str(pv_accuracy)]);
disp(['The pv PEFAC has precision: ', num2str(pv_precision)]);
disp(['The pv PEFAC has recall: ', num2str(pv_recall)]);

disp([' ']);
disp([' ']);

% print the metrics
disp(['The new pv method has accuracy: ', num2str(x_m_accuracy)]);
disp(['The new pv method has precision: ', num2str(x_m_precision)]);
disp(['The new pv method has recall: ', num2str(x_m_recall)]);

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Things to do: 
% 1) EM algorithm VS LS/MMSE for training the KF. We use only LS/MMSE here. We need to use the EM algorithm.

% 2) Use future frames. Use KF for smoothing/hindsight and not for filtering. In general, filtering VS prediction VS smoothing.

% 3) Now, pv is the frame voiced probability. We need to try pv(k,l), we need to try fr.bin and frame voiced probability.

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% 1)
% KF training: LS way VS EM algorithm

% previously, we defined total_total_observations 
% previously, we defined the states: total_states

% we now need to train the KF again but with the EM algorithm
% KF -- train parameters -- assume constant

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% 2)
% we use fixed-lag smoothing

% we use future frames 
% we now use the forward-backward algorithm

% we use noisy speech signal
% we implement the KF filter

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

reference_true_pv = pv;

% Read the noise file
[s2, fs2] = readwav('babble.wav');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

predicted_x_m = [];
predicted_P_m = [];

% we define the lag d
d = 4;

for loop = 1 : length(pv)
    % we initialize x_m and delta_delta_phi and P_p
    if (loop == 1)
        x_m = 0.5;
        P_m = 1000.2;
    end
    
    % define l 
    % l is the observation
    l = pv(loop);

    K = ((F*P_m*F') + Q) * H' * ((H * ((F*P_m*F') + Q) * H') + R)^(-1);

    x_p = (F*x_m) + (K * (l - (H * F * x_m)));

    P_p = (1-K) * ((F * P_m * F') + Q);

    if (loop > d && (loop+d) <= length(pv))
        future_x_p = [];
        future_P_p = [];
        
        for kk1 = 1 : d
            % define l 
            % l is the observation
            l = pv(loop+kk1);

            K = ((F*P_m*F') + Q) * H' * ((H * ((F*P_m*F') + Q) * H') + R)^(-1);

            x_p = (F*x_m) + (K * (l - (H * F * x_m)));

            P_p = (1-K) * ((F * P_m * F') + Q);
            
            future_x_p = [future_x_p x_p];
            future_P_p = [future_P_p P_p]; 
        end

        % we can use weighted sum
        % or we can use the mean value
        
        %x_p = mean([x_p future_x_p]);
        x_p = (0.5*x_p) + (0.4*future_x_p(1));
        
        for kk2 = 2 : d
            x_p = x_p + ((0.1/(d-2)) * future_x_p(kk2));
        end
        
        %P_p = mean([P_p future_P_p]);
        P_p = (0.5*P_p) + (0.4*future_P_p(1));
        
        for kk2 = 2 : d
            P_p = P_p + ((0.1/(d-2)) * future_P_p(kk2));
        end
        
    end
    
    % constrain x_p
    if (x_p < 0)
       x_p = 0.01; 
    end

    % constrain x_p
    if (x_p > 1)
       x_p = 0.99; 
    end

    % store the results
    predicted_x_m = [predicted_x_m x_p];
    predicted_P_m = [predicted_P_m P_p];

    % prepare the next iteration
    x_m = x_p;
    P_m = P_p;

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);

plot(predicted_x_m);
hold on;
plot(pv, 'g');
hold on;
plot(reference_true_pv, 'r');
hold off;

set(gca,'FontSize',10);
title('pv. Probability of frame being voiced.');
xlabel('Time in frames');
ylabel('Probability pv');

legend(' new pv', ' PEFAC pv', ' true pv');

% ------------------------------------------------------------------------

% treat the problem as a classification problem
% we use 0 and 1 only

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) > 0.5)
       reference_true_pv(i) = 1;
   else
       reference_true_pv(i) = 0;
   end
end

for i = 1 : length(pv)
   if (pv(i) > 0.5)
       pv(i) = 1;
   else
       pv(i) = 0;
   end
end

for i = 1 : length(predicted_x_m)
   if (predicted_x_m(i) > 0.5)
       predicted_x_m(i) = 1;
   else
       predicted_x_m(i) = 0;
   end
end

% we find accuracy for pv
pv_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == pv(i))
       pv_accuracy = pv_accuracy + 1;
   end
end

pv_accuracy = pv_accuracy / length(reference_true_pv);

% we find accuracy for x_m
x_m_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == predicted_x_m(i))
       x_m_accuracy = x_m_accuracy + 1;
   end
end

x_m_accuracy = x_m_accuracy / length(reference_true_pv);

% we now use pv
% we compare pv and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

pv_true_positives = 0;
pv_true_negatives = 0;

pv_false_positives = 0;
pv_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && pv(i) == 1)
        pv_true_positives = pv_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 0)
        pv_true_negatives = pv_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && pv(i) == 0)
        pv_false_positives = pv_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 1)
        pv_false_negatives = pv_false_negatives + 1;
    end
end
    
% find precision
pv_precision = (pv_true_positives) / (pv_true_positives+pv_false_positives);

% find recall
pv_recall = (pv_true_positives) / (pv_true_positives+pv_false_negatives);

% we now do the same for x_m

% we now use x_m
% we compare x_m and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

x_m_true_positives = 0;
x_m_true_negatives = 0;

x_m_false_positives = 0;
x_m_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && predicted_x_m(i) == 1)
        x_m_true_positives = x_m_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 0)
        x_m_true_negatives = x_m_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && predicted_x_m(i) == 0)
        x_m_false_positives = x_m_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 1)
        x_m_false_negatives = x_m_false_negatives + 1;
    end
end

% find precision
x_m_precision = (x_m_true_positives) / (x_m_true_positives+x_m_false_positives);

% find recall
x_m_recall = (x_m_true_positives) / (x_m_true_positives+x_m_false_negatives);

% we thus have 3 metrics for the two methods

% we have: pv_accuracy, pv_precision, pv_recall 
% we have: x_m_accuracy, x_m_precision, x_m_recall

disp([' ']);
disp([' ']);

disp([' ']);
disp([' ']);

% print the metrics
disp(['The pv PEFAC has accuracy: ', num2str(pv_accuracy)]);
disp(['The pv PEFAC has precision: ', num2str(pv_precision)]);
disp(['The pv PEFAC has recall: ', num2str(pv_recall)]);

disp([' ']);
disp([' ']);

% print the metrics
disp(['The new pv method has accuracy: ', num2str(x_m_accuracy)]);
disp(['The new pv method has precision: ', num2str(x_m_precision)]);
disp(['The new pv method has recall: ', num2str(x_m_recall)]);

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


adfsadfasfa

















%%
% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');
 
% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);
 
% use PEFAC
[total_harmonics_total_freq, total_harmonics_total_amplitude] = activlevg_new_new(s,fs);
 
% we have 200 voiced frames
 
% we have in time: total_observations
% we have harmonics frequency: total_harmonics_total_freq
% we have harmonics amplitude: total_harmonics_total_amplitude
 
% we use total_harmonics_total_freq
new_main_total_harmonics_total_freq = total_harmonics_total_freq;
new_main_total_harmonics_total_freq(end+1,:) = new_main_total_harmonics_total_freq(end,:);
 
% we use pv
voiced = [];
for i = 1 : length(pv)
    if(pv(i) > 0.5)
        voiced(i) = 1;
    else
        voiced(i) = 0;
    end
end
 
% Define the frame size in samples
% len = txinc_samples*10;
len = 0.090 * fs;
 
% window overlap in percent of frame size
PERC = (80/90) * 100;
 
len1 = floor(len * PERC/100);
len2 = len - len1;
 
segment_shift_L = len2;
 
x = s;
 
% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;
 
nFFT = 2 * len;
 
k = 1;
img = sqrt(-1);
 
% Define window
win = hamming(len);
 
% we use NFFT
disp('This is the total frequency bins:')
nFFT 
 
% we use unique NFFT
disp('This is the unique frequency bins:')
nFFT_unique = (nFFT/2) + 1;
nFFT_unique
 
% find the omega_h
omega_h = (total_harmonics_total_freq/fs) * 2 * pi;
 
% find omega_k
k = 0:1:nFFT_unique-1;
omega_k = 2 * pi* (k) * (1/nFFT_unique);
 
% initialize the omega_k_h
omega_k_h = zeros(size(omega_h,1), length(omega_k));
omega_k_h2 = zeros(size(omega_h,1), length(omega_k));

% we set the function that we will minimize
for outter_loop = 1 : size(omega_h,1)
    for inner_loop = 1 : length(omega_k)
        array_temp_value_minimum = [];
        
        clear k;
        
        for k = 1 : size(omega_h,2)
            array_temp_value_minimum = [array_temp_value_minimum abs(omega_k(inner_loop)-omega_h(outter_loop,k))];
        end
        
        value_index = find(array_temp_value_minimum==min(array_temp_value_minimum));
        value_index = value_index(1);
        
        temp_value_minimum = omega_h(outter_loop, value_index);
        
        % we want to set the omega_h value correctly
        % we want to find the minimum value
        omega_k_h(outter_loop, inner_loop) = temp_value_minimum;
        
        % find the index
        
        % find the f in omega_k
        for kk = 1 : length(omega_k)
           if (omega_k(kk) <= temp_value_minimum && omega_k(kk+1) >= temp_value_minimum)
               main_index_value = kk;
           end
        end
        
        omega_k_h2(outter_loop, inner_loop) = main_index_value;
    
    end
end
 
% we now have omega_k_h
 
final_omega_k_h = omega_k_h;

% we use omega_k_h
omega_k_h_new = zeros(size(omega_k_h));
omega_k_h(end+1,:) = omega_k_h(end,:);

omega_k_h2_new = zeros(size(omega_k_h2));
omega_k_h2(end+1,:) = omega_k_h2(end,:);

for pp = 1 : length(voiced)
    if(voiced(pp) == 1)
         omega_k_h_new(pp,:) = omega_k_h(1,:);
         omega_k_h = omega_k_h(2:end,:);
         
         omega_k_h2_new(pp,:) = omega_k_h2(1,:);
         omega_k_h2 = omega_k_h2(2:end,:);

    else
        omega_k_h_new(pp,:) = zeros(size(omega_k_h(1,:)));
        
        omega_k_h2_new(pp,:) = zeros(size(omega_k_h2(1,:)));

    end
end

% double-sided
double_sided_omega_k_h = [flipud(omega_k_h_new(:,3:end)) omega_k_h_new];
 
double_sided_omega_k_h2 = [flipud(omega_k_h2_new(:,3:end)) omega_k_h2_new];

% initialize variable
first_time_voiced = 1;
 
% phase_estimate_clean_speech = [];
phase_estimate_clean_speech = zeros(Nframes,nFFT);
 
% define N_W
N_w_window = 620;

omega_k_h = final_omega_k_h;

f_h_k2 = omega_k_h * fs * (1/(2*pi));
f_h_k2(end+1,:) = f_h_k2(end,:);

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    signal_part = x(k:k+len-1);
    insign = signal_part .* win;
 
    % Take Fourier transform of frame
    spec = fft(insign, nFFT);
 
    insign_phase = angle(spec);
    
    one_sided_insign_phase = insign_phase(((end/2)+1):end);
    
    % Take Fourier transform of window
    spec2 = fft(win, nFFT);
 
    insign_phase2 = angle(spec2);

    one_sided_insign_phase2 = insign_phase2(((end/2)+1):end);
    
    % we use 3 cases
    % cases: 1) first voiced, 2) not first voiced, 3) non-voiced frames
    
    % ------------------------------------------------------------------------
    % ------------------------------------------------------------------------

    % this is case 1
    % when we find the first voiced frame
    
    if (voiced(n) == 1 && first_time_voiced == 1)
        first_time_voiced = 0;
        
        f_h_k = f_h_k2(1,:);
        f_h_k2 = f_h_k2(2:end,:);

        %temp_phase_estimate_clean_speech = one_sided_insign_phase';
        temp_phase_estimate_clean_speech = insign_phase';
        
        % we only keep the fr.Bins with harmonics
        % we disregard the other fr.Bins
        
        % ------------------------------------------------------------------------
        % ------------------------------------------------------------------------

        % when fr.Bins that contain harmonics
        
        % we use new_main_total_harmonics_total_freq(1,:)
        current_new_main_total_harmonics_total_freq = new_main_total_harmonics_total_freq(1,:);
        new_main_total_harmonics_total_freq = new_main_total_harmonics_total_freq(2:end,:);
        
        omega_current_new_main_total_harmonics_total_freq = (2*pi) * (1/fs) * current_new_main_total_harmonics_total_freq;
        
        % delete the zero elements
        omega_current_new_main_total_harmonics_total_freq(omega_current_new_main_total_harmonics_total_freq==0) = [];
        
%         % one-sided fr.Bins
%         omega_one_sided = (2*pi) * (0:1:(nFFT_unique-1)) * (1/fs);
%         
%         % two-sided fr.Bins
%         omega_two_sided = [flipud(omega_one_sided(2:end)) omega_one_sided];
        
        % we use double_sided_omega_k_h
        for loop_var1 = 1 : nFFT
            %harmonic_to_be_used = a(n,loop_var1);
            
            % omega_k = [flipud(omega_k(2:end)) omega_k];
            
            %index1 = find(omega_k == harmonic_to_be_used);
            %index1 = index1(1);
            index1 = double_sided_omega_k_h2(n,loop_var1);
            
            phase_estimate_clean_speech(n,loop_var1) = temp_phase_estimate_clean_speech(index1);
        
        end
        
        % we only keep the fr.Bins with harmonics
        % we disregard the other fr.Bins
        
        for loop_var1 = 2 : nFFT
            if (phase_estimate_clean_speech(n,loop_var1) == phase_estimate_clean_speech(n,loop_var1-1))
               phase_estimate_clean_speech(n,loop_var1) = 0;
               
            end
        end
            
        % ------------------------------------------------------------------------
        % ------------------------------------------------------------------------

        % when fr.Bins that do not contain harmonics
        
%         for loop_var1 = 1 : nFFT
%             harmonic_to_be_used = double_sided_omega_k_h(n,loop_var1);
%             
%             %index1 = find(omega_k == harmonic_to_be_used);
%             index1 = double_sided_omega_k_h2(n,loop_var1);
%             
%             if (phase_estimate_clean_speech(n,loop_var1) == 0)
%                 % we need to find window_dependent_var1
% 
%                 % we use insign_phase2
%                 window_dependent_var1 = ((-1) * insign_phase2(round(index1-(N_w_window*(1/fs)*f_h_k(loop_var1))))) + (insign_phase2(index1+(loop_var1-index1)-(N_w_window*(1/fs)*f_h_k(loop_var1))));
% 
%                 phase_estimate_clean_speech(n,loop_var1) = phase_estimate_clean_speech(n,index1) + (window_dependent_var1);
%             
%             end 
%         end

        for loop_var1 = (nFFT/2)+1 : nFFT
            harmonic_to_be_used = double_sided_omega_k_h(n,loop_var1);
            
            %index1 = find(omega_k == harmonic_to_be_used);
            index1 = double_sided_omega_k_h2(n,loop_var1);
            
            if (phase_estimate_clean_speech(n,loop_var1) == 0)
                % we need to find window_dependent_var1

                % we use insign_phase2
                window_dependent_var1 = ((-1) * insign_phase2(abs(round(index1-(N_w_window*(1/fs)*f_h_k(loop_var1-(nFFT/2))))))) + (insign_phase2(round(index1+abs(loop_var1-index1)-(N_w_window*(1/fs)*f_h_k(loop_var1-(nFFT/2))))));

                phase_estimate_clean_speech(n,loop_var1) = phase_estimate_clean_speech(n,index1) + (window_dependent_var1);
                phase_estimate_clean_speech(n,1+(nFFT/2)-(loop_var1-(nFFT/2))) = phase_estimate_clean_speech(n,loop_var1);
            
            end 
        end

    % ------------------------------------------------------------------------
    % ------------------------------------------------------------------------
    
    % this is case 2
    % not first voiced frame
    
    elseif (voiced(n) == 1 && first_time_voiced == 0)
        f_h_k2 = f_h_k2(2:end,:);

        % we use f_h_k
        % we use f_h_k(loop1)
        
        for loop1 = (nFFT/2)+1 : nFFT
            phase_estimate_clean_speech(n,loop1) = phase_estimate_clean_speech(n-1,loop1) + ((2*pi*len2)*(1/fs)*(f_h_k(loop1-(nFFT/2))));
            phase_estimate_clean_speech(n,1+(nFFT/2)-(loop_var1-(nFFT/2))) = phase_estimate_clean_speech(n,loop_var1);
            
        end
        
    % ------------------------------------------------------------------------
    % ------------------------------------------------------------------------
    
    % this is case 3
    % non-voiced frame
    
    elseif (voiced(n) == 0)
        first_time_voiced = 1;
        phase_estimate_clean_speech(n,:) = zeros(size(insign_phase'));
        
    end
 
    k = k+len2;
end

% we found phase_estimate_clean_speech
% we estimated the STFT clean phase

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


afsafsafdas













%%
% Kulmer 
% we implement Kulmer's paper for phase estimation

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% Read the noise file
[s2, fs2, bits] = readwav('babble.wav');

% Define the SNR
snr = 15;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we use pv
voiced = [];
for i = 1 : length(pv)
    if(pv(i) > 0.5)
        voiced(i) = 1;
    else
        voiced(i) = 0;
    end
end

% use PEFAC
[total_harmonics_total_freq, total_harmonics_total_amplitude] = activlevg_new_new_final(s,fs);

% frame increment [fs=sample frequency]
ninc = 0.010 * fs;   

% overlap factor
% PERC = 100 * (80/90);
% ovf=1/PERC;                  
ovf = 1440/ninc;

f=rfft(enframe(s,hamming(ovf*ninc,'periodic'),ninc),2*ovf*ninc,2);

% convert to power spectrum
f=f.*conj(f);           

% estimate the noise power spectrum
noise_estimate_frBins_frame = estnoisem(f,ninc/fs); 

final_total_harmonics_total_freq = [];
final_total_harmonics_total_amplitude = [];

for loop = 1 : length(voiced)
    if (voiced(loop) == 1)
        final_total_harmonics_total_freq = [final_total_harmonics_total_freq; total_harmonics_total_freq(1,:)];
        final_total_harmonics_total_amplitude = [final_total_harmonics_total_amplitude; total_harmonics_total_amplitude(1,:)];
        
        total_harmonics_total_freq = total_harmonics_total_freq(2:end,:);
        total_harmonics_total_amplitude = total_harmonics_total_amplitude(2:end,:);
        
    else
        final_total_harmonics_total_freq = [final_total_harmonics_total_freq; zeros(1,size(total_harmonics_total_amplitude,2))];
        final_total_harmonics_total_amplitude = [final_total_harmonics_total_amplitude; zeros(1,size(total_harmonics_total_amplitude,2))];
        
    end
end

% we use final_total_harmonics_total_freq
% we use final_total_harmonics_total_amplitude

% Define the frame size in samples
% len = txinc_samples*10;
len = 0.090 * fs;

% window overlap in percent of frame size
PERC = (80/90) * 100;

len1 = floor(len * PERC/100);
len2 = len - len1;

segment_shift_L = len2;

x = s;

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;

nFFT = 2 * len;

% we use NFFT
disp('This is the total frequency bins:')
nFFT 

% we use unique NFFT
disp('This is the unique frequency bins:')
nFFT_unique = (nFFT/2) + 1;
nFFT_unique

k = 1;
img = sqrt(-1);

% Define window
win = hamming(len);

% % plot win
% figure; 
% set(gcf,'Color','w');
% set(gca,'FontSize',10);
% plot(win);
% set(gca,'FontSize',10);
% title('Window.');
% xlabel('Samples, Time.');
% ylabel('Amplitude.');

N_w_window = 620;

N = len;

% find the omega
final_total_harmonics_total_freq = (final_total_harmonics_total_freq/fs) * 2 * pi;

pv = (pv/fs) * 2 * pi;

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    signal_part = x(k:k+len-1);
    insign = signal_part .* win;

    % we only use the voiced frames

    % we use noise estimate
    % we use noise_estimate_frBins_frame
    
    % we use fx from PEFAC
    % fx is the pitch for every frame
    
    % detect the voiced frames
    if (voiced(n) == 1)
        % we use final_total_harmonics_total_freq
        % we use final_total_harmonics_total_amplitude
        
        current_frame_final_total_harmonics_total_freq = final_total_harmonics_total_freq(n,:);
        current_frame_final_total_harmonics_total_amplitude = final_total_harmonics_total_amplitude(n,:);

        current_frame_noise_estimate_frBins_frame = noise_estimate_frBins_frame(n,:);

        % we need to convert the noise estimate from frBins to harmonics

        % define omeaga_h
        omega_h = current_frame_final_total_harmonics_total_freq;
        
        % find omega_k
        k = 0:1:nFFT_unique-1;
        omega_k = 2 * pi* (k) * (1/nFFT_unique);

        % initialize the omega_k_h
        omega_k_h = [];
        omega_k_h = zeros(size(omega_h,1), length(omega_k));

        % we set the function that we will minimize
        for outter_loop = 1 : size(omega_h,1)
            for inner_loop = 1 : length(omega_k)
                array_temp_value_minimum = [];

                clear k;

                for k = 1 : size(omega_h,2)
                    array_temp_value_minimum = [array_temp_value_minimum abs(omega_k(inner_loop)-omega_h(outter_loop,k))];
                end

                value_index = find(array_temp_value_minimum==min(array_temp_value_minimum));
                value_index = value_index(1);

                temp_value_minimum = omega_h(outter_loop, value_index);

                % we want to set the omega_h value correctly
                % we want to find the minimum value
                omega_k_h(outter_loop, inner_loop) = temp_value_minimum;

            end
        end

        % we now have omega_k_h
        
        % we use omega_k_h and current_frame_noise_estimate_frBins_frame
        current_frame_noise_estimate_h_frame = current_frame_noise_estimate_frBins_frame(1);
        
        for ttt = 2 : length(current_frame_noise_estimate_frBins_frame)
            if (omega_k_h(ttt) ~= omega_k_h(ttt-1))
                current_frame_noise_estimate_h_frame = [current_frame_noise_estimate_h_frame current_frame_noise_estimate_frBins_frame(ttt)];
                
            end
        end
                
        current_frame_noise_estimate_h_frame_final = [];
        
        for u = 1 : length(current_frame_final_total_harmonics_total_freq)
            if (current_frame_final_total_harmonics_total_freq(u) == 0)
                current_frame_noise_estimate_h_frame_final = [current_frame_noise_estimate_h_frame_final 0];
            
            else
                current_frame_noise_estimate_h_frame_final = [current_frame_noise_estimate_h_frame_final current_frame_noise_estimate_h_frame(1)];
                current_frame_noise_estimate_h_frame = current_frame_noise_estimate_h_frame(2:end);
                
            end
        end
                
        % we found the noise
        noise_h_l = current_frame_noise_estimate_h_frame_final;     
        
        % the phase phi_h_l needs to be found

        h = 0 : 1 : 14;
        sum1 = zeros(size(h));

        for jj = 1 : N
           sum1 = sum1 + (insign(jj) * sin(h*fx(n)*jj)); 
        end
        
        von_mises = 0;
        numerator1 = (-((2*current_frame_final_total_harmonics_total_amplitude)./noise_h_l).*sum1) + (von_mises);
        
        h = 0 : 1 : 14;
        sum2 = zeros(size(h));

        for jj = 1 : N
           sum2 = sum2 + (insign(jj) * cos(h*fx(n)*jj)); 
        end
        
        von_mises2 = 0.2;
        denominator1 = (-((2*current_frame_final_total_harmonics_total_amplitude)./noise_h_l).*sum2) + (von_mises2);

        phi_h_l = atan(numerator1 ./ denominator1);
        
        % the phase phi_h_l is found
        % we need the phase phi_k,l
        
%         % phi_h_l and phi_k_l are different
%         
%         N_p = min(N_w_window, ((fx(n)*nFFT_unique)/(2*pi)));
%         
%         % i = -(N_p/2) : 1 : (N_p/2);
%                 
%         % fx is the pitch for every frame
%         h_harmonic = 0:1:14;
%         A_variable = floor(h_harmonic * fx(n) * nFFT_unique) + i;
%         
%         phi_k_l = phi_h_l(A_variable);
        
        phi_k_l = zeros(size(omega_k_h));
        phi_k_l(1) = phi_h_l(1);
        
        for pp = 2 : length(omega_k_h)
            if (omega_k_h(pp) == omega_k_h(pp-1))
                phi_k_l(pp) = phi_k_l(pp-1);
                
            else
                phi_k_l(pp) = phi_h_l(1);
                phi_h_l = phi_h_l(2:end);
                
            end
            
        end
    end

    k = k+len2;
end

% we implemented Kulmer's paper for phase estimation

% we found phi_k_l
% we used phi_k_l and phi_h_l 

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

asfasdfasfdasfda







%%
% we fit Wrapped GMMs (WGMMs) to phase histogram
% we fit von-Mises distribution to phase histogram

% we use the Kalman Filter (KF) for clean phase estimation

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we use pv
voiced = [];
for i = 1 : length(pv)
    if(pv(i) > 0.5)
        voiced(i) = 1;
    else
        voiced(i) = 0;
    end
end

% we only keep voiced frames

% Define the frame size in samples
% len = txinc_samples*10;
len = 0.090 * fs;

% window overlap in percent of frame size
PERC = (80/90) * 100;

len1 = floor(len * PERC/100);
len2 = len - len1;

segment_shift_L = len2;

x = s;

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;

nFFT = 2 * len;

k = 1;
img = sqrt(-1);

% Define window
win = hamming(len);

total_observations = [];
total_observations_phase = [];
total_observations_phase_one_side = [];

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    signal_part = x(k:k+len-1);
    insign = signal_part .* win;

    % Take Fourier transform of  frame
    spec = fft(insign, nFFT);

    insign_phase = angle(spec);
    
    if (voiced(n) == 1)
        total_observations = [total_observations; insign'];
        
        total_observations_phase_one_side = [total_observations_phase_one_side; insign_phase(length(insign_phase)/2:end)'];
        
        total_observations_phase = [total_observations_phase; insign_phase'];
    end

    k = k+len2;
end

% we have 200 voiced frames
% we have in time: total_observations

% we have phase: total_observations_phase

% we choose one frequency bin
% we model all frequency bins with one frequency bin

total_observations_phase_one_bin = total_observations_phase(:, floor(size(total_observations_phase,2)/2));
total_observations_phase_one_bin_one_side = total_observations_phase_one_side(:, floor(size(total_observations_phase,2)/2));

% this was the first signal
% we now try the second signal

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA2.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we use pv
voiced = [];
for i = 1 : length(pv)
    if(pv(i) > 0.5)
        voiced(i) = 1;
    else
        voiced(i) = 0;
    end
end

% we only keep voiced frames

% Define the frame size in samples
% len = txinc_samples*10;
len = 0.090 * fs;

% window overlap in percent of frame size
PERC = (80/90) * 100;

len1 = floor(len * PERC/100);
len2 = len - len1;

segment_shift_L = len2;

x = s;

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;

nFFT = 2 * len;

k = 1;
img = sqrt(-1);

% Define window
win = hamming(len);

total_observations = [];
total_observations_phase = [];
total_observations_phase_one_side = [];

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    signal_part = x(k:k+len-1);
    insign = signal_part .* win;

    % Take Fourier transform of  frame
    spec = fft(insign, nFFT);

    insign_phase = angle(spec);
    
    if (voiced(n) == 1)
        total_observations = [total_observations; insign'];
        
        total_observations_phase_one_side = [total_observations_phase_one_side; insign_phase(length(insign_phase)/2:end)'];
        
        total_observations_phase = [total_observations_phase; insign_phase'];
    end

    k = k+len2;
end

% we have 200 voiced frames
% we have in time: total_observations

% we have phase: total_observations_phase

% we choose one frequency bin
% we model all frequency bins with one frequency bin

total_observations_phase_one_bin2 = total_observations_phase(:, floor(size(total_observations_phase,2)/2));

total_observations_phase_one_bin = [total_observations_phase_one_bin; total_observations_phase_one_bin2];

total_observations_phase_one_bin_one_side = [total_observations_phase_one_bin_one_side; total_observations_phase_one_side(:, floor(size(total_observations_phase,2)/2))];

[s,fs,wrd,phn] = readsph('SI648.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we use pv
voiced = [];
for i = 1 : length(pv)
    if(pv(i) > 0.5)
        voiced(i) = 1;
    else
        voiced(i) = 0;
    end
end

% we only keep voiced frames

% Define the frame size in samples
% len = txinc_samples*10;
len = 0.090 * fs;

% window overlap in percent of frame size
PERC = (80/90) * 100;

len1 = floor(len * PERC/100);
len2 = len - len1;

segment_shift_L = len2;

x = s;

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;

nFFT = 2 * len;

k = 1;
img = sqrt(-1);

% Define window
win = hamming(len);

total_observations = [];
total_observations_phase = [];
total_observations_phase_one_side = [];

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    signal_part = x(k:k+len-1);
    insign = signal_part .* win;

    % Take Fourier transform of  frame
    spec = fft(insign, nFFT);

    insign_phase = angle(spec);
    
    if (voiced(n) == 1)
        total_observations = [total_observations; insign'];
        
        total_observations_phase_one_side = [total_observations_phase_one_side; insign_phase(length(insign_phase)/2:end)'];
        
        total_observations_phase = [total_observations_phase; insign_phase'];
    end

    k = k+len2;
end

% we have 200 voiced frames
% we have in time: total_observations

% we have phase: total_observations_phase

% we choose one frequency bin
% we model all frequency bins with one frequency bin

total_observations_phase_one_bin3 = total_observations_phase(:, floor(size(total_observations_phase,2)/2));

total_observations_phase_one_bin = [total_observations_phase_one_bin; total_observations_phase_one_bin3];

total_observations_phase_one_bin_one_side = [total_observations_phase_one_bin_one_side; total_observations_phase_one_side(:, floor(size(total_observations_phase,2)/2))];

% % plot the phase histogram 
% figure; 
% set(gcf,'Color','w');
% set(gca,'FontSize',10);
% subplot(2,2,1);
% hist(total_observations_phase_one_bin, 50);
% set(gca,'FontSize',10);
% title('Histogram.');
% xlabel('Phase');
% ylabel('Count');
% 
% % plot the log phase histogram 
% figure; 
% set(gcf,'Color','w');
% set(gca,'FontSize',10);
% subplot(2,2,1);
% hist(log(total_observations_phase_one_bin+5), 50);
% set(gca,'FontSize',10);
% title('Histogram.');
% xlabel('Log Phase. log(phase+5).');
% ylabel('Count');

% phase is between 0 and 2*pi 
% or phase is between -pi and pi

% we use: phase between -pi and pi

for i = 1 : length(total_observations_phase_one_bin)
    if (total_observations_phase_one_bin(i) > pi)
        total_observations_phase_one_bin(i) = total_observations_phase_one_bin(i) - (2*pi);
        
    elseif (total_observations_phase_one_bin(i) < -pi)
        total_observations_phase_one_bin(i) = total_observations_phase_one_bin(i) + (2*pi);
                
    end

end

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin, 50);
set(gca,'FontSize',10);
title('Histogram. Phase between -pi and pi.');
xlabel('Phase');
ylabel('Count');

% plot the log phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(log(total_observations_phase_one_bin+5), 50);
set(gca,'FontSize',10);
title('Histogram. Phase between -pi and pi.');
xlabel('Log Phase. log(phase+5).');
ylabel('Count');

% we now use the one sided phase

% phase is between 0 and 2*pi 
% or phase is between -pi and pi

% we use: phase between -pi and pi

for i = 1 : length(total_observations_phase_one_bin_one_side)
    if (total_observations_phase_one_bin_one_side(i) > pi)
        total_observations_phase_one_bin_one_side(i) = total_observations_phase_one_bin_one_side(i) - (2*pi);
        
    elseif (total_observations_phase_one_bin_one_side(i) < -pi)
        total_observations_phase_one_bin_one_side(i) = total_observations_phase_one_bin_one_side(i) + (2*pi);
                
    end

end

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
set(gca,'FontSize',10);
title('Histogram. One Side Phase between -pi and pi.');
xlabel('Phase');
ylabel('Count');

% plot the log phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(log(total_observations_phase_one_bin_one_side+5), 50);
set(gca,'FontSize',10);
title('Histogram. One Side Phase between -pi and pi.');
xlabel('Log Phase. log(phase+5).');
ylabel('Count');

% for -pi to pi, we use the von Mises pd
phase_phi_s = -pi : 0.1 : pi;

kappa_k = 0.2;

von_mises_pd = (1/(2*pi)) * (1/besseli(0,kappa_k)) * exp(kappa_k * cos(phase_phi_s));

figure;
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
plot(phase_phi_s, von_mises_pd);
set(gca,'FontSize',10);
title('Von Mises pd.');
xlabel('RV');
ylabel('pd');

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
hold on;
plot(phase_phi_s, von_mises_pd*100, 'r');
hold off;

set(gca,'FontSize',10);
title('Histogram. One Side Phase between -pi and pi.');
xlabel('Phase');
ylabel('Count');

% for -pi to pi, we use the von Mises pd
phase_phi_s = -pi : 0.1 : pi;
kappa_k = 4.82;
von_mises_pd = (1/(2*pi)) * (1/besseli(0,kappa_k)) * exp(kappa_k * cos(phase_phi_s));

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
hold on;
plot(phase_phi_s, von_mises_pd*50, 'r');
hold off;

set(gca,'FontSize',10);
title('Histogram. Kappa = 4.82. One Side Phase between -pi and pi.');
xlabel('Phase');
ylabel('Count');

% we now use: phase between 0 and 2*pi

for i = 1 : length(total_observations_phase_one_bin)
    if (total_observations_phase_one_bin(i) > 2*pi)
        total_observations_phase_one_bin(i) = total_observations_phase_one_bin(i) - (2*pi);
        
    elseif (total_observations_phase_one_bin(i) < 0)
        total_observations_phase_one_bin(i) = total_observations_phase_one_bin(i) + (2*pi);
                
    end
end

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin, 50);
set(gca,'FontSize',10);
title('Histogram. Phase between 0 and 2pi.');
xlabel('Phase');
ylabel('Count');

% plot the log phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(log(total_observations_phase_one_bin+5), 50);
set(gca,'FontSize',10);
title('Histogram. Phase between 0 and 2pi.');
xlabel('Log Phase. log(phase+5).');
ylabel('Count');

% we use the one sided

% we now use: phase between 0 and 2*pi

for i = 1 : length(total_observations_phase_one_bin_one_side)
    if (total_observations_phase_one_bin_one_side(i) > 2*pi)
        total_observations_phase_one_bin_one_side(i) = total_observations_phase_one_bin_one_side(i) - (2*pi);
        
    elseif (total_observations_phase_one_bin_one_side(i) < 0)
        total_observations_phase_one_bin_one_side(i) = total_observations_phase_one_bin_one_side(i) + (2*pi);
                
    end
end

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
set(gca,'FontSize',10);
title('Histogram. One Side Phase between 0 and 2pi.');
xlabel('Phase');
ylabel('Count');

% plot the log phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(log(total_observations_phase_one_bin_one_side+5), 50);
set(gca,'FontSize',10);
title('Histogram. One Side Phase between 0 and 2pi.');
xlabel('Log Phase. log(phase+5).');
ylabel('Count');

% when 0 to 2*pi, we use the wrapped Gaussian
% we use the wrapped Gaussian pd

x = 0 : 0.1 : 2*pi;
f_x_prob = zeros(size(x));

mean_mu = pi;
variance_sigma_squared = 0.2;

for k = -1 : 1 : 1
    f_x_prob = f_x_prob + ((1/sqrt(2*pi*variance_sigma_squared)) * (exp(-((x+(2*pi*k)-mean_mu).^2)/(2*variance_sigma_squared))));
    
end

% plot the wrapped Gaussian pd

figure;
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
plot(x, f_x_prob);
set(gca,'FontSize',10);
title('Wrapped Gaussian pd. mu = pi and variance = 0.2.');
xlabel('RV');
ylabel('pd');

% we now use the EM algorithm

% we can use EM to fit the wrapped Gaussian pd
% we use training set for the EM algorithm

% we use initial mu_m, sigma_squared_s values

% % we use random initializations
% mu_m = pi;
% sigma_squared_s = 0.5;

% total_observations_phase_one_bin_one_side = total_observations_phase_one_bin_one_side(randperm(length(total_observations_phase_one_bin_one_side)));

[m,v] = gaussmix(total_observations_phase_one_bin_one_side,[],[],1);

mu_m = m;
sigma_squared_s = v;

% we do 10 iterations
for iterations_EM = 1 : 20

    % we use P_x_k
    % k is from -1 to 1 

    % we use total_observations_phase_one_bin_one_side

    total_number_of_examples = length(total_observations_phase_one_bin_one_side);

    P_x_k_total_total = [];

    for x_i = 1 : total_number_of_examples
        x_current_example = total_observations_phase_one_bin_one_side(x_i);

        P_x_k_total = [];
        
        % do the E step
        for k = -1 : 1 : 1
            denominator = 0;
            x = x_current_example;

            for k = -1 : 1 : 1
                denominator = denominator + ((1/sqrt(2*pi*sigma_squared_s)) * (exp(-((x+(2*pi*k)-mean_mu).^2)/(2*sigma_squared_s))));

            end

            % to avoid dividing by zero
            denominator = denominator + 0.001;
            
            P_x_k = ((1/(sqrt(2*pi*sigma_squared_s))) * (exp(-(x_current_example+(2*pi*k-mu_m))/(2*sigma_squared_s)))) / (denominator);

            P_x_k_total = [P_x_k_total P_x_k];

        end
        
        P_x_k_total_total = [P_x_k_total_total; P_x_k_total];
    end

    % we use sum(P_x_k_total_total')
    % we use mean(sum(P_x_k_total_total'))
    
    matrix_temp = [(total_observations_phase_one_bin_one_side-(2*pi)) (total_observations_phase_one_bin_one_side) (total_observations_phase_one_bin_one_side+(2*pi))];
    
    P_x_k_total_total2 = P_x_k_total_total .* matrix_temp;
    
    mu_m = mean(sum(P_x_k_total_total2'));
    
    mu_m = mod(mu_m,(2*pi));
    
    % we use sum(P_x_k_total_total')
    % we use mean(sum(P_x_k_total_total'))
    
    matrix_temp = [(total_observations_phase_one_bin_one_side-(2*pi)-mu_m).^2 (total_observations_phase_one_bin_one_side-mu_m).^2 (total_observations_phase_one_bin_one_side+(2*pi)-mu_m).^2];
    
    P_x_k_total_total2 = P_x_k_total_total .* matrix_temp;
    
    sigma_squared_s = mean(sum(P_x_k_total_total2'));
end

% we now have the wrapped Gaussian pd

% we plot the wrapped Gaussian pd

x = 0 : 0.1 : 2*pi;
f_x_prob = zeros(size(x));

mean_mu = mu_m;
variance_sigma_squared = sigma_squared_s;

for k = -1 : 1 : 1
    f_x_prob = f_x_prob + ((1/sqrt(2*pi*variance_sigma_squared)) * (exp(-((x+(2*pi*k)-mean_mu).^2)/(2*variance_sigma_squared))));
    
end

% plot the wrapped Gaussian pd

figure;
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
plot(x, f_x_prob);
set(gca,'FontSize',10);
title('Wrapped Gaussian pd. mu and variance were estimated.');
xlabel('RV');
ylabel('pd');

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
hold on;
plot(x, f_x_prob*50, 'r');
hold off;

set(gca,'FontSize',10);
title('Histogram. mu and variance were estimated. One Side Phase between 0 and 2*pi.');
xlabel('Phase');
ylabel('Count');

% % we use M Gaussians
% M = 2;
% 
% psi_n = total_observations_phase_one_bin_one_side;
% 
% delta_m_n_minus_2
% delta_m_n_minus_1
% delta_m_n_zero
% delta_m_n_plus_1
% delta_m_n_plus_2
% 
% beta_m_n = 0;
% 
% for w = -2 : 1 : 2
%     beta_m_n = beta_m_n + ();
%     
% end
% 
% beta_m_n_mu = 0;
% 
% for w = -2 : 1 : 2
%     beta_m_n_mu = beta_m_n_mu + ();
%     
% end
% 
% beta_m_n_sigma_squared = 0;
% 
% for w = -2 : 1 : 2
%     beta_m_n_sigma_squared = beta_m_n_sigma_squared + ();
%     
% end
% 
% beta_m_n_new = a_m * beta_m_n;
% 
% w_m = 0;
% N = length(total_observations_phase_one_bin_one_side);
% 
% for n = 1 : N
%     temp_var1 = 0;
%     
%     for u = 1 : M
%         temp_var1 = temp_var1 + (beta_m_n);
%     end
%     
%     w_m = w_m + (beta_m_n/temp_var1);    
%     
% end
% 
% a_m = (1/N) * w_m;
% 
% sum1 = 0;
% 
% for n = 1 : N
%     sum1 = sum1 + ((beta_m_n_new/beta_m_n) * beta_m_n_mu);
% end
% 
% m = (1/w_m) * (sum1);

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% we repeat the above process
% we now use k from -2 to 2

% we now use the EM algorithm

% we can use EM to fit the wrapped Gaussian pd
% we use training set for the EM algorithm

% we use initial mu_m, sigma_squared_s values

% % we use random initializations
% mu_m = pi;
% sigma_squared_s = 0.5;

% total_observations_phase_one_bin_one_side = total_observations_phase_one_bin_one_side(randperm(length(total_observations_phase_one_bin_one_side)));

[m,v] = gaussmix(total_observations_phase_one_bin_one_side,[],[],1);

mu_m = m;
sigma_squared_s = v;

% we do 10 iterations
for iterations_EM = 1 : 20

    % we use P_x_k
    % k is from -1 to 1 

    % we use total_observations_phase_one_bin_one_side

    total_number_of_examples = length(total_observations_phase_one_bin_one_side);

    P_x_k_total_total = [];

    for x_i = 1 : total_number_of_examples
        x_current_example = total_observations_phase_one_bin_one_side(x_i);

        P_x_k_total = [];
        
        % do the E step
        for k = -2 : 1 : 2
            denominator = 0;
            x = x_current_example;

            for k = -2 : 1 : 2
                denominator = denominator + ((1/sqrt(2*pi*sigma_squared_s)) * (exp(-((x+(2*pi*k)-mean_mu).^2)/(2*sigma_squared_s))));

            end

            % to avoid dividing by zero
            denominator = denominator + 0.001;
            
            P_x_k = ((1/(sqrt(2*pi*sigma_squared_s))) * (exp(-(x_current_example+(2*pi*k-mu_m))/(2*sigma_squared_s)))) / (denominator);

            P_x_k_total = [P_x_k_total P_x_k];

        end
        
        P_x_k_total_total = [P_x_k_total_total; P_x_k_total];
    end

    % we use sum(P_x_k_total_total')
    % we use mean(sum(P_x_k_total_total'))
    
    matrix_temp = [(total_observations_phase_one_bin_one_side-(4*pi)) (total_observations_phase_one_bin_one_side-(2*pi)) (total_observations_phase_one_bin_one_side) (total_observations_phase_one_bin_one_side+(2*pi)) (total_observations_phase_one_bin_one_side+(4*pi))];
    
    P_x_k_total_total2 = P_x_k_total_total .* matrix_temp;
    
    mu_m = mean(sum(P_x_k_total_total2'));
    
    % we use sum(P_x_k_total_total')
    % we use mean(sum(P_x_k_total_total'))
    
    matrix_temp = [(total_observations_phase_one_bin_one_side-(4*pi)-mu_m).^2 (total_observations_phase_one_bin_one_side-(2*pi)-mu_m).^2 (total_observations_phase_one_bin_one_side-mu_m).^2 (total_observations_phase_one_bin_one_side+(2*pi)-mu_m).^2 (total_observations_phase_one_bin_one_side+(4*pi)-mu_m).^2];
    
    P_x_k_total_total2 = P_x_k_total_total .* matrix_temp;
    
    sigma_squared_s = mean(sum(P_x_k_total_total2'));
end

% we now have the wrapped Gaussian pd

% we plot the wrapped Gaussian pd

x = 0 : 0.1 : 2*pi;
f_x_prob = zeros(size(x));

mean_mu = mu_m;
variance_sigma_squared = sigma_squared_s;

for k = -2 : 1 : 2
    f_x_prob = f_x_prob + ((1/sqrt(2*pi*variance_sigma_squared)) * (exp(-((x+(2*pi*k)-mean_mu).^2)/(2*variance_sigma_squared))));
    
end

% plot the wrapped Gaussian pd

figure;
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
plot(x, f_x_prob);
set(gca,'FontSize',10);
title('Wrapped Gaussian pd. mu and variance were estimated.');
xlabel('RV');
ylabel('pd');

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
hold on;
plot(x, f_x_prob*50, 'r');
hold off;

set(gca,'FontSize',10);
title('Histogram. mu and variance were estimated. One Side Phase between 0 and 2*pi.');
xlabel('Phase');
ylabel('Count');

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now use another way

% we use alternative way
theta_n = total_observations_phase_one_bin_one_side;

mean_circular = angle(mean(exp(sqrt(-1)*theta_n)));
sigma_squared_circular = 1 - abs(mean(exp(sqrt(-1)*theta_n)));

final_mean_wrapped_Gaussian = mod(mean_circular, 2*pi);
final_sigma_squared_wrapped_Gaussian = -2 * log(1-sigma_squared_circular);

% we now have the wrapped Gaussian pd

% we plot the wrapped Gaussian pd

x = 0 : 0.1 : 2*pi;
f_x_prob = zeros(size(x));

mean_mu = final_mean_wrapped_Gaussian;
variance_sigma_squared = final_sigma_squared_wrapped_Gaussian;

for k = -1 : 1 : 1
    f_x_prob = f_x_prob + ((1/sqrt(2*pi*variance_sigma_squared)) * (exp(-((x+(2*pi*k)-mean_mu).^2)/(2*variance_sigma_squared))));
    
end

% plot the wrapped Gaussian pd

figure;
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
plot(x, f_x_prob);
set(gca,'FontSize',10);
title('Wrapped Gaussian pd. mu and variance were estimated.');
xlabel('RV');
ylabel('pd');

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
hold on;
plot(x, f_x_prob*50, 'r');
hold off;

set(gca,'FontSize',10);
title('Histogram. mu and variance were estimated. One Side Phase between 0 and 2*pi.');
xlabel('Phase');
ylabel('Count');

% % we investigate the relation of gaussmix
% theta_n = total_observations_phase_one_bin_one_side;
% 
% [m,v] = gaussmix(theta_n,[],[],1);
% 
% % we use the wrapped Gaussian pd
% 
% x = 0 : 0.1 : 2*pi;
% f_x_prob = zeros(size(x));
% 
% constant1 = (6.2016/3.1809);
% mean_mu = m*constant1;
% 
% constant2 = 1.3421/5.8370;
% variance_sigma_squared = v*constant2;
% 
% for k = -1 : 1 : 1
%     f_x_prob = f_x_prob + ((1/sqrt(2*pi*variance_sigma_squared)) * (exp(-((x+(2*pi*k)-mean_mu).^2)/(2*variance_sigma_squared))));
%     
% end
% 
% % plot the wrapped Gaussian pd
% 
% figure;
% set(gcf,'Color','w');
% set(gca,'FontSize',10);
% subplot(2,2,1);
% plot(x, f_x_prob);
% set(gca,'FontSize',10);
% title('Wrapped Gaussian pd. We use gaussmix for mu and variance.');
% xlabel('RV');
% ylabel('pd');

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now use another way

f_total = [];
final_prob_total = [];

for loop = 1 : 10
    M = 2;
    [m,v,w,~,f]=gaussmix(total_observations_phase_one_bin_one_side,[],[],M,'kf');

    weights_w = w;

    mu = mod(m, (2*pi));

    prob = 0 : 0.1 : (2*pi);

    final_prob = zeros(size(prob));

    for m = 1 : M
        sum1 = 0;

        for w = -2 : 1 : 2
            sum1 = sum1 + ((1/(sqrt(2*pi*v(m)))) * (exp(-(1/(2*v(m)))*(prob-mu(m)-(2*pi*w)).^2)));

        end

        final_prob = final_prob + (weights_w(m) * sum1);    

    end
    
    f_total = [f_total f];
    final_prob_total = [final_prob_total; final_prob];
end

% find the largest f
[~, index1] = max(f_total);

% f_total = f_total(index1);
final_prob_total = final_prob_total(index1,:);

% plot the wrapped Gaussian pd

figure;
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
plot(prob, final_prob);
set(gca,'FontSize',10);
title('Wrapped Gaussian pd. We use gaussmix for mu and variance.');
xlabel('RV');
ylabel('pd');

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
hold on;
plot(prob, final_prob*50, 'r');
hold off;

set(gca,'FontSize',10);
title('Histogram. mu and variance were estimated. One Side Phase between 0 and 2*pi.');
xlabel('Phase');
ylabel('Count');

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now use 3 components in the GMM

f_total = [];
final_prob_total = [];

for loop = 1 : 10
    M = 3;
    [m,v,w,~,f]=gaussmix(total_observations_phase_one_bin_one_side,[],[],M,'kf');

    weights_w = w;

    mu = mod(m, (2*pi));

    prob = 0 : 0.1 : (2*pi);

    final_prob = zeros(size(prob));

    for m = 1 : M
        sum1 = 0;

        for w = -2 : 1 : 2
            sum1 = sum1 + ((1/(sqrt(2*pi*v(m)))) * (exp(-(1/(2*v(m)))*(prob-mu(m)-(2*pi*w)).^2)));

        end

        final_prob = final_prob + (weights_w(m) * sum1);    

    end
    
    f_total = [f_total f];
    final_prob_total = [final_prob_total; final_prob];
end

% find the largest f
[~, index1] = max(f_total);

% f_total = f_total(index1);
final_prob_total = final_prob_total(index1,:);

% plot the phase histogram 
figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);
subplot(2,2,1);
hist(total_observations_phase_one_bin_one_side, 50);
hold on;
plot(prob, final_prob*50, 'r');
hold off;

set(gca,'FontSize',10);
title('Histogram. mu and variance were estimated. One Side Phase between 0 and 2*pi.');
xlabel('Phase');
ylabel('Count');

% % -----------------------------------------------------------------------
% % -----------------------------------------------------------------------
% % -----------------------------------------------------------------------
% 
% % let us now consider the KF
% 
% % KF for an AR model
% x_n = randn(1,1);
% 
% a_1 = 0.78;
% 
% for i = 1 : 999
%    x_n = [x_n (randn(1,1)-(a_1*x_n(end)))]; 
% end
% 
% % plot the AR(1) process
% 
% figure; 
% set(gcf,'Color','w');
% set(gca,'FontSize',10);
% subplot(2,2,1);
% plot(1:length(x_n), x_n);
% 
% set(gca,'FontSize',10);
% title('AR(1) process.');
% xlabel('Time');
% ylabel('Amplitude');

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% training the KF

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we use pv
voiced = [];
for i = 1 : length(pv)
    if(pv(i) > 0.5)
        voiced(i) = 1;
    else
        voiced(i) = 0;
    end
end

% we only keep voiced frames

% Define the frame size in samples
% len = txinc_samples*10;
len = 0.090 * fs;

% window overlap in percent of frame size
PERC = (80/90) * 100;

len1 = floor(len * PERC/100);
len2 = len - len1;

segment_shift_L = len2;

x = s;

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;

nFFT = 2 * len;

k = 1;
img = sqrt(-1);

% Define window
win = hamming(len);

total_observations = [];
real_STFT_phases = [];

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    signal_part = x(k:k+len-1);
    insign = signal_part .* win;

    % Take Fourier transform of  frame
    spec = fft(insign, nFFT);

    insign_phase = angle(spec);
    
    if (voiced(n) == 1)
        total_observations = [total_observations; insign'];
        
        real_STFT_phases = [real_STFT_phases; insign_phase'];
    
    else
        real_STFT_phases = [real_STFT_phases; zeros(size(insign_phase))'];
    
    end

    k = k+len2;
end

% use PEFAC
total_harmonics_total_freq = activlevg_new_new(s,fs);

% we have 200 voiced frames

% we have in time: total_observations
% we have harmonics frequency: total_harmonics_total_freq
% we have harmonics amplitude: total_harmonics_total_amplitude

% we need observations 
% we need pitch differences

pitch_difference = total_harmonics_total_freq(:,1);

pitch_difference2 = zeros(size(voiced));

for j_loop = 1 : length(pitch_difference2)
    if (voiced(j_loop) == 1)
        pitch_difference2(j_loop) = pitch_difference(1);
        
        pitch_difference = pitch_difference(2:end);
    
    else
        pitch_difference2(j_loop) = 0;
        
    end
end

pitch_difference = pitch_difference2;

pitch_difference = diff(pitch_difference);

total_total_observations = pitch_difference';

% we need states
% we need real STFT phases

% we use real_STFT_phases
real_STFT_phases = real_STFT_phases';

% KF -- train parameters -- assume constant

% In KF, there is no access to the correct data when the algorithm is running in real time.
% In KF, the error is between the measurement and the Kalman prediction.
% It can be stated that the KF resembles a low pass filter.
% but its transfer characteristic is nonlinear, and the cut-off frequency shifts.

% In KF, there are 2 main things to do in each time step. 
% The first thing is to calculate the apriori state estimate, and the measurement update.
% The second step is to calculate the gain, corresponding covariance, and then the posterior state.

total_states = real_STFT_phases;
total_states = total_states(:,2:end);

total_states = total_states((size(total_states,1)/2+22),:);

% we use M time steps
M = size(total_states,2);

for loop = 1 : M
    temp = total_states(:,loop);
    
    X(:,loop) = temp';
end

X1 = X(:, 1:(end-1));
X2 = X(:, 2:end);

clear temp;
C = 1;

total_total_observations = total_total_observations';

for loop = 1 : M
        temp = zeros(1,C);

        % temp = observations;
        temp = total_total_observations(:,loop);

        Z(:,loop) = temp';
end

%A = X2 * X1' * (X1 * X1')^(-1);
A = X2 * X1' * pinv(X1 * X1');

%H = Z * X' * (X * X')^(-1);
H = Z * X' * pinv(X * X');

W = (X2 - (A * X1)) * (X2 - (A * X1))' / (M-1);

Q = (Z - (H * X)) * (Z - (H * X))' / M; 

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% the true phases are real_STFT_phases

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% do the KF

% A is the same as the F matrix
% W is the same as the Q matrix

F = A;
Q1 = Q;
Q = W;

% H is the same as H
% Q is the same as the R matrix

R = Q1;

% we use noisy speech signal
% we implement the KF filter

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% Read the noise file
[s2, fs2, bits] = readwav('babble.wav');

% Define the SNR
snr = 15;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we use pv
voiced = [];
for i = 1 : length(pv)
    if(pv(i) > 0.5)
        voiced(i) = 1;
    else
        voiced(i) = 0;
    end
end

% we only keep voiced frames

% Define the frame size in samples
% len = txinc_samples*10;
len = 0.090 * fs;

% window overlap in percent of frame size
PERC = (80/90) * 100;

len1 = floor(len * PERC/100);
len2 = len - len1;

segment_shift_L = len2;

x = s;

% Define the parameters
Nframes = floor(length(x)/len2)-1;
Nframes = Nframes - 7;

nFFT = 2 * len;

k = 1;
img = sqrt(-1);

% Define window
win = hamming(len);

first_time = 1;

estimated_phases = [];

% Main for loop
for n = 1 : Nframes
    % Define the input signal
    signal_part = x(k:k+len-1);
    insign = signal_part .* win;

    % Take Fourier transform of  frame
    spec = fft(insign, nFFT);

    insign_phase = angle(spec);

    if (n == 1)
        phase_estimate = zeros(size(insign_phase));
    end
    
    if (voiced(n) == 1)
        % we do not use KF
        if (first_time == 1)
            phase_estimate = insign_phase;
            
            first_time = 0;
        else
            % we now use KF
            observation_y = fx(n) - fx(n-1);

            x_m = phase_estimate;

            x_derivative_m = F * x_m;
            x_m = x_derivative_m;

            % initialize delta_delta_phi
            delta_delta_phi = rand(size(Q));

            delta_delta_t_phi = F * delta_delta_phi;
            delta_delta_phi = delta_delta_t_phi;

            % initialize P_p
            P_p = rand(size(Q))*10;

            P_m = (delta_delta_phi * P_p * (delta_delta_phi)') + Q;
            K = P_m * H' * (H*P_m*H'+R)^(-1);
            l_m = H * x_m;

            % define l 
            l = observation_y;

            x_p = x_m + (K * (l-l_m));

            I_matrix = eye(size(K*H));
            P_p = (I_matrix - (K*H)) * P_m;

            % prepare the next iteration
            x_m = x_p;
            phase_estimate = x_m;
            P_m = P_p;
        end
        
        estimated_phases = [estimated_phases; phase_estimate'];
        
    else
        first_time = 1;
        
        estimated_phases = [estimated_phases; zeros(size(phase_estimate))'];
    end

    k = k + len2;
end

% compare with the true phases
% we use real_STFT_phases 

% we use estimated_phases

percentage_error = mean(abs(estimated_phases(:)-real_STFT_phases(:)) ./ (real_STFT_phases(:)+0.001));

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% we do not do the KF with the wrapped GMM

% we now have to try the wrapped GMM

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
