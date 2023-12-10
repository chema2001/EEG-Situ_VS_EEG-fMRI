%% Falta fazer:

% - Análises estatística dos picos 
% - Ajustar os thesholds dos picos e ver melhor literatura dos picos entre o
% N170/VPP e P300
% Se der tempo, aplicar ICA

%% Informação dos picos

% N170 (100-200 ms) - Potencial com amplitude maior durante a apresenção de faces em
% comparação com outros objetos/estímulos visuais. Elétrodos: P7, P8, PO7, PO8.
% Geralmente é mais pronunciando no hemisfério direito.
% Vertex Positive Potential (VPP) - 100 - 200 ms. Parece ser gerado pelo mesmo 
% processo neuronal de reconhecimento facil mas é captado nos elétrodos centrais: 
% FZ, CZ, apesar do poder haver pequenas diferenças de latência.
% Mais sobre o N170 e VPP: https://psyc.bbk.ac.uk/brainb/wp-content/uploads/sites/13/2019/05/eimer_N170_oxford.pdf

% N2 (200-400 ms): FZ, Cz
% P3 (300-600 ms): FZ, Fz, Pz
% https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1469-8986.2007.00602.x

%% Importar os dados e inicializar variáveis

close all;
clc, clear;

fig = 1;

% EEGLAB functions documentation: https://sccn.ucsd.edu/~arno/eeglab/auto/indexfunc.html

addpath(genpath('eeglab2023.0')); % Path do EEGLAB 
ch_locs=readlocs('BioSemi64.loc'); % Ficheiro com localização dos elétrodos

data_seated = readtable("S3_25_Male.csv");
data_seated = table2array(data_seated);

data_seatedsound = readtable("S2_25_sentado_som_2023-06-21_18-18.csv");
data_seatedsound = table2array(data_seatedsound);

data_fmri = readtable("S2_25_Male_2023-06-02_16-21.csv");
data_fmri = table2array(data_fmri);

data_lying = readtable("S22_25_Male_2023-06-02_17-03.csv");
data_lying = table2array(data_lying);

fs = 125;
channel_labels = ["C3";"CP3";"P3";"P7";"PO7";"PO3";"Pz";"CPz";"Fz";"Cz";"C4";"CP4";"P4";"P8";"PO8";"PO4"];

idx = zeros(length(channel_labels),1);

for i = 1:length(ch_locs)
    cur_lb = string(ch_locs(i).labels);
    if (nonzeros(channel_labels == cur_lb))
        idx(i) = 1;
    end
end

nch_locs = ch_locs(logical(idx)); % Elétrodos usados

%% 2D plot


%figure (fig), topoplot([], nch_locs, 'style', 'blank', 'electrodes', 'ptslabels');
%title(["2D Channel Locations";""]);

%fig = fig + 1;

%% Headplot

splfile='FileHeadPlot.spl';
%headplot('setup',nch_locs,splfile,'Param','Value');
headplot(data_fmri(100000,1:end-1),splfile,'electrodes','on','labels',2,'maplimits', 'absmax','colormap', sky(64));

%% Verificar o efeito dos filtros na fase dos sinais

% t = (-1:1/fs:1)';
% s = 5*square(2*pi*t); % Sinal com a mesma frequência que o sinal da fotoresistência
% 
% sqwave = awgn(s, 20, 'measured'); % Gaussian white noise
% 
% plot(t, [sqwave, movmean(sqwave,20), movmedian(sqwave,20)]),
% title('Filter effects on the phase of the phoresistor signal'), 
% xlabel('Time (s)'),
% legend('Original Signal','Moving average','Moving median');

% Resumindo, a mediana flutuante não distorce tanto as fases do sinal,
% principalmente no início das rising edges e falling edges. Isto é
% extremamente importante em análises como as dos ERPs pois é preciso
% saber com exatidão o instante em que o estímulo é apresentado.

%% Dados da fotoresistência 

photoVector_seated=(data_seated(:,end));

t_seated = ((0:length(photoVector_seated)-1)*1/fs)/60;
n_pnts = 20;
filt_photoVector_seated = movmedian(photoVector_seated,20);

photoVector_seatedsound=(data_seatedsound(:,end));

t_seatedsound = ((0:length(photoVector_seatedsound)-1)*1/fs)/60;
filt_photoVector_seatedsound = movmedian(photoVector_seatedsound,20);

% lpFilt = designfilt('lowpassfir','FilterOrder',1000,'PassbandFrequency',2, ...
%          'StopbandFrequency',2.5,'SampleRate',fs);
%fvtool(lpFilt)

%lpFilt_phVec = filtfilt(lpFilt,photoVector);

photoVector_fmri=(data_fmri(:,end));

t_fmri = ((0:length(photoVector_fmri)-1)*1/fs)/60;
filt_photoVector_fmri = movmedian(photoVector_fmri,20);

photoVector_lying=(data_lying(:,end));

t_lying = ((0:length(photoVector_lying)-1)*1/fs)/60;
filt_photoVector_lying = movmedian(photoVector_lying,20);

%figure (fig); plot(t,photoVector), title('Fotoresistência'), xlabel('Tempo (min)'); % Sem filtros
%fig = fig + 1;

figure (fig); plot(t_seated,filt_photoVector_seated), title("Fotoresistência (Mediana flutuante com "+n_pnts+" pontos) - EEG Sentado"), xlabel('Tempo (min)');
fig = fig + 1;

figure (fig); plot(t_seatedsound,filt_photoVector_seatedsound), title("Fotoresistência (Mediana flutuante com "+n_pnts+" pontos) - EEG Sentado com som"), xlabel('Tempo (min)');
fig = fig + 1;

% figure (fig); plot(t,lpFilt_phVec), title("Fotoresistência (Filtro passa-baixo)"), xlabel('Tempo (min)');
% fig = fig + 1;

figure (fig); plot(t_fmri,filt_photoVector_fmri), title("Fotoresistência (Mediana flutuante com "+n_pnts+" pontos) - EEG com fMRI"), xlabel('Tempo (min)');
fig = fig + 1;

figure (fig); plot(t_lying,filt_photoVector_lying), title("Fotoresistência (Mediana flutuante com "+n_pnts+" pontos) - EEG Deitado"), xlabel('Tempo (min)');
fig = fig + 1;

%% Re-referenciação do EEG 

eeg_data_seated = data_seated(:,1:end-1);
eeg_data_seatedsound = data_seatedsound(:,1:end-1);
eeg_data_fmri = data_fmri(:,1:end-1);
eeg_data_lying = data_lying(:,1:end-1);

%car = bsxfun(@minus, eeg_data, mean(eeg_data,2));
%lap = laplacian_perrinX(eeg_data',[nch_locs.X],[nch_locs.Y],[nch_locs.Z]);

%% Pré-processamento do EEG

% Filtragem
% De acordo com o artigo com a literatura dos ERPs, um
% filtro passa-banda entre 0.5 e 30/40 Hz será usado.

% Band-pass
bpFilt = designfilt('bandpassfir','FilterOrder',3000, ...
         'CutoffFrequency1',0.1,'CutoffFrequency2',30, ...
         'SampleRate',fs); 
%fvtool(bpFilt) % Visualizar resposta em frequência e fase do filtro

%filt_data = filtfilt(bpFilt,car);
%filt_data = filtfilt(bpFilt,lap');

filt_data_seated = filtfilt(bpFilt,eeg_data_seated); % Aplica ZeroPhase
filt_data_seatedsound = filtfilt(bpFilt,eeg_data_seatedsound);
filt_data_fmri = filtfilt(bpFilt,eeg_data_fmri);
filt_data_lying = filtfilt(bpFilt,eeg_data_lying);

%% FFT do sinal 

% Verificar se o sinal é um EEG. Para isso, o sinal no domíno das
% frequências tem que seguir a forma 1/f.
% Isto após retirar as componentes DC e dos 50Hz da rede elétrica.

% ----- Sentado -----

T = 1/fs; % Sampling period
L = length(filt_data_seated(:,12)); % Length of signal
tm = (0:L-1)*T; % Time vector
Y = fft(filt_data_seated(:,8)); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:floor(L/2)+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L; % Frequency vector

figure(fig)
fig = fig +1;
plot(f,P1) 
title('Espectro de Amplitudes X(t) - Sentado')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% ----- Sentado som -----

T = 1/fs; % Sampling period
L = length(filt_data_seatedsound(:,8)); % Length of signal
tm = (0:L-1)*T; % Time vector
Y = fft(filt_data_seatedsound(:,8)); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:floor(L/2)+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L; % Frequency vector

figure(fig)
fig = fig +1;
plot(f,P1) 
title('Espectro de Amplitudes X(t) - Sentado com som')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% ----- fMRI -----

T = 1/fs; % Sampling period
L = length(filt_data_fmri(:,8)); % Length of signal
tm = (0:L-1)*T; % Time vector
Y = fft(filt_data_fmri(:,8)); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:floor(L/2)+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L; % Frequency vector

figure(fig)
fig = fig +1;
plot(f,P1) 
title('Espectro de Amplitudes X(t) - EEG com fMRI')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% ----- Deitado -----

T = 1/fs; % Sampling period
L = length(filt_data_lying(:,8)); % Length of signal
tm = (0:L-1)*T; % Time vector
Y = fft(filt_data_lying(:,8)); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:floor(L/2)+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L; % Frequency vector

figure(fig)
fig = fig +1;
plot(f,P1) 
title('Espectro de Amplitudes X(t) - Deitado')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Retirar sub-harmónica dos 25 Hz

%Band-stop
bsFilt = designfilt('bandstopfir','FilterOrder',3000, ...
         'CutoffFrequency1',24.8,'CutoffFrequency2',25.2, ...
         'SampleRate',fs);
%fvtool(bsFilt)

filt_data_seated = filtfilt(bsFilt,filt_data_seated);
filt_data_seatedsound = filtfilt(bsFilt,filt_data_seatedsound);
filt_data_fmri = filtfilt(bsFilt,filt_data_fmri);
filt_data_lying = filtfilt(bsFilt,filt_data_lying);

% ICA

%filt_data_seated = pop_loadset()


%% Deteção dos índices de segmentação

% ----- Seated -----

last_low = true; 
last_high = false;

falling_edges = [];
rising_edges = [];

for i = 1:length(filt_photoVector_seated)
    photoresistor = filt_photoVector_seated(i);
    if photoresistor <= 235 % Possible falling edge
        if last_high % Found falling edge
            falling_edges = [falling_edges;i];
        end

        last_low = true; 
        last_high = false;

    elseif photoresistor >= 235 % Possible rising edge
        if last_low % Found rising edge
            rising_edges = [rising_edges; i];
        end
        last_low = false;
        last_high = true; 
    end
end

rising_edges(1) = []; 

seg_target_seated = [];
seg_nontarget_seated = [];

for i = 1:length(falling_edges) % Verificar quais blocos pertecem aos targets e não targets
    
    aux = filt_photoVector_seated(falling_edges(i):rising_edges(i));

    if (length(aux)>=60) && (length(aux)<=85)
    
        if(min(aux) < 140)
            seg_nontarget_seated = [seg_nontarget_seated, falling_edges(i)];
        else
            seg_target_seated = [seg_target_seated, falling_edges(i)];
        end

    end

end

seg_target_seated = unique(seg_target_seated);
seg_nontarget_seated = unique(seg_nontarget_seated);


% ----- Seated w/ sound -----

last_low = true; 
last_high = false;

falling_edges = [];
rising_edges = [];

for i = 1:length(filt_photoVector_seatedsound)
    photoresistor = filt_photoVector_seatedsound(i);
    if photoresistor <= 232 % Possible falling edge
        if last_high % Found falling edge
            falling_edges = [falling_edges;i];
        end

        last_low = true; 
        last_high = false;

    elseif photoresistor >= 232 % Possible rising edge
        if last_low % Found rising edge
            rising_edges = [rising_edges; i];
        end
        last_low = false;
        last_high = true; 
    end
end

falling_edges(1:2) = [];
rising_edges(1:3) = [];

falling_edges(end) = [];

seg_target_seatedsound = [];
seg_nontarget_seatedsound = [];

for i = 1:length(falling_edges) % Verificar quais blocos pertecem aos targets e não targets
    
    aux = filt_photoVector_seatedsound(falling_edges(i):rising_edges(i));

    if (length(aux)>=50) && (length(aux)<=75)
    
        if(min(aux) < 130)
            seg_nontarget_seatedsound = [seg_nontarget_seatedsound, falling_edges(i)];
        else
            seg_target_seatedsound = [seg_target_seatedsound, falling_edges(i)];
        end

    end

end

seg_target_seatedsound = unique(seg_target_seatedsound);
seg_nontarget_seatedsound = unique(seg_nontarget_seatedsound);

% ----- fMRI -----

last_low = true; 
last_high = false;

falling_edges = [];
rising_edges = [];

for i = 1:length(filt_photoVector_fmri)
    photoresistor = filt_photoVector_fmri(i);
    if photoresistor <= 233 % Possible falling edge
        if last_high % Found falling edge
            falling_edges = [falling_edges;i];
        end

        last_low = true; 
        last_high = false;

    elseif photoresistor >= 233 % Possible rising edge
        if last_low % Found rising edge
            rising_edges = [rising_edges; i];
        end
        last_low = false;
        last_high = true; 
    end
end

rising_edges(1) = [];
falling_edges(end) = [];

seg_target_fmri = [];
seg_nontarget_fmri = [];

for i = 1:length(falling_edges) % Verificar quais blocos pertecem aos targets e não targets
    
    aux = filt_photoVector_fmri(falling_edges(i):rising_edges(i));
    
    if (length(aux)>=60) && (length(aux)<=80)

        if(min(aux) < 130)
            seg_nontarget_fmri = [seg_nontarget_fmri, falling_edges(i)];
        else
            seg_target_fmri = [seg_target_fmri, falling_edges(i)];
        end

    end
end

seg_target_fmri = unique(seg_target_fmri);
seg_nontarget_fmri = unique(seg_nontarget_fmri);

% ----- Deitado -----

last_low = true; 
last_high = false;

falling_edges = [];
rising_edges = [];

for i = 1:length(filt_photoVector_lying)
    photoresistor = filt_photoVector_lying(i);
    if photoresistor <= 235 % Possible falling edge
        if last_high % Found falling edge
            falling_edges = [falling_edges;i];
        end

        last_low = true; 
        last_high = false;

    elseif photoresistor >= 230 % Possible rising edge
        if last_low % Found rising edge
            rising_edges = [rising_edges; i];
        end
        last_low = false;
        last_high = true; 
    end
end

rising_edges(1) = [];
falling_edges(end) = [];

seg_target_lying = [];
seg_nontarget_lying = [];

for i = 1:length(falling_edges) % Verificar quais blocos pertecem aos targets e não targets
    
    aux = filt_photoVector_lying(falling_edges(i):rising_edges(i));
    
    if (length(aux)>=60) && (length(aux)<=80)

        if(min(filt_photoVector_lying(falling_edges(i):rising_edges(i))) < 130)
            seg_nontarget_lying = [seg_nontarget_lying, falling_edges(i)];
        else
            seg_target_lying = [seg_target_lying, falling_edges(i)];
        end

    end
end

seg_target_lying = unique(seg_target_lying);
seg_nontarget_lying = unique(seg_nontarget_lying);

%% Segmentação do EEG

bf = -0.200; % 100 ms antes do estímulo
af = 0.800; % 800 ms após o estímulo

idx_bf = floor(bf * fs);
idx_af = ceil(af * fs);

time = ((idx_bf:idx_af)/fs) * 1000; % Vetor do tempo

% -------------------
% Target - Sentado

[targetEEG_seated, nontargetEEG_seated] = deal(zeros(length(seg_target_seated), length(time), length(channel_labels))); % Pré-alocação dos dados

for i=1:length(seg_target_seated)
    targetEEG_seated(i,:,:) = cat(3,filt_data_seated((seg_target_seated(i)+idx_bf):(seg_target_seated(i)+idx_af),:));
end

% Non-target - Sentado

for i=1:length(seg_nontarget_seated)
    nontargetEEG_seated(i,:,:) = cat(3,filt_data_seated((seg_nontarget_seated(i)+idx_bf):(seg_nontarget_seated(i)+idx_af),:));
end

% -------------------
% Target - Sentado c/ som

[targetEEG_seatedsound, nontargetEEG_seatedsound] = deal(zeros(length(seg_target_seatedsound), length(time), length(channel_labels))); % Pré-alocação dos dados

for i=1:length(seg_target_seatedsound)
    targetEEG_seatedsound(i,:,:) = cat(3,filt_data_seatedsound((seg_target_seatedsound(i)+idx_bf):(seg_target_seatedsound(i)+idx_af),:));
end

% Non-target - Sentado c/ som

for i=1:length(seg_nontarget_seatedsound)
    nontargetEEG_seatedsound(i,:,:) = cat(3,filt_data_seatedsound((seg_nontarget_seatedsound(i)+idx_bf):(seg_nontarget_seatedsound(i)+idx_af),:));
end

% -------------------
% Target - fMRI

[targetEEG_fmri, nontargetEEG_fmri] = deal(zeros(length(seg_target_fmri), length(time), length(channel_labels))); % Pré-alocação dos dados

for i=1:length(seg_target_fmri)
    targetEEG_fmri(i,:,:) = cat(3,filt_data_fmri((seg_target_fmri(i)+idx_bf):(seg_target_fmri(i)+idx_af),:));
end

% Non-target - fMRI

for i=1:length(seg_nontarget_fmri)
    nontargetEEG_fmri(i,:,:) = cat(3,filt_data_fmri((seg_nontarget_fmri(i)+idx_bf):(seg_nontarget_fmri(i)+idx_af),:));
end

% -------------------
% Target - Deitado

[targetEEG_lying, nontargetEEG_lying] = deal(zeros(length(seg_target_lying), length(time), length(channel_labels))); % Pré-alocação dos dados

for i=1:length(seg_target_lying)
    targetEEG_lying(i,:,:) = cat(3,filt_data_lying((seg_target_lying(i)+idx_bf):(seg_target_lying(i)+idx_af),:));
end

% Non-target - Deitado

for i=1:length(seg_nontarget_lying)
    nontargetEEG_lying(i,:,:) = cat(3,filt_data_lying((seg_nontarget_lying(i)+idx_bf):(seg_nontarget_lying(i)+idx_af),:));
end


%% Correção da linha de base

baseline = [-200;0];

baseidx = dsearchn(time',baseline);

% ----- Sentado -----

for t=1:size(targetEEG_seated,1)
    for e=1:size(targetEEG_seated,3)
        tmn = mean(targetEEG_seated(t,baseidx(1):baseidx(2),e));   
        targetEEG_seated(t,:,e) = targetEEG_seated(t,:,e) - tmn;
    end
end 

for t=1:size(nontargetEEG_seated,1)
    for e=1:size(nontargetEEG_seated,3)
        tmn = mean(nontargetEEG_seated(t,baseidx(1):baseidx(2),e));   
        nontargetEEG_seated(t,:,e) = nontargetEEG_seated(t,:,e) - tmn;
    end
end 

% ----- Sentado c/som -----

for t=1:size(targetEEG_seatedsound,1)
    for e=1:size(targetEEG_seatedsound,3)
        tmn = mean(targetEEG_seatedsound(t,baseidx(1):baseidx(2),e));   
        targetEEG_seatedsound(t,:,e) = targetEEG_seatedsound(t,:,e) - tmn;
    end
end 

for t=1:size(nontargetEEG_seatedsound,1)
    for e=1:size(nontargetEEG_seatedsound,3)
        tmn = mean(nontargetEEG_seatedsound(t,baseidx(1):baseidx(2),e));   
        nontargetEEG_seatedsound(t,:,e) = nontargetEEG_seatedsound(t,:,e) - tmn;
    end
end 

% ----- fMRI -----

for t=1:size(targetEEG_fmri,1)
    for e=1:size(targetEEG_fmri,3)
        tmn = mean(targetEEG_fmri(t,baseidx(1):baseidx(2),e));   
        targetEEG_fmri(t,:,e) = targetEEG_fmri(t,:,e) - tmn;
    end
end 

for t=1:size(nontargetEEG_fmri,1)
    for e=1:size(nontargetEEG_fmri,3)
        tmn = mean(nontargetEEG_fmri(t,baseidx(1):baseidx(2),e));   
        nontargetEEG_fmri(t,:,e) = nontargetEEG_fmri(t,:,e) - tmn;
    end
end 

% ----- Deitado -----

for t=1:size(targetEEG_lying,1)
    for e=1:size(targetEEG_lying,3)
        tmn = mean(targetEEG_lying(t,baseidx(1):baseidx(2),e));   
        targetEEG_lying(t,:,e) = targetEEG_lying(t,:,e) - tmn;
    end
end 

for t=1:size(nontargetEEG_lying,1)
    for e=1:size(nontargetEEG_lying,3)
        tmn = mean(nontargetEEG_lying(t,baseidx(1):baseidx(2),e));   
        nontargetEEG_lying(t,:,e) = nontargetEEG_lying(t,:,e) - tmn;
    end
end 

%% Computação dos ERPs

% ----- Sentado -----

tERP_seated = squeeze(mean(targetEEG_seated, 1)); % ERP target (média do sinal na dimensão trials)
ntERP_seated = squeeze(mean(nontargetEEG_seated, 1)); % ERP não-target

% ----- Sentado c/ som -----

tERP_seatedsound = squeeze(mean(targetEEG_seatedsound, 1)); % ERP target (média do sinal na dimensão trials)
ntERP_seatedsound = squeeze(mean(nontargetEEG_seatedsound, 1)); % ERP não-target


% ----- fMRI -----

tERP_fmri = squeeze(mean(targetEEG_fmri, 1)); % ERP target (média do sinal na dimensão trials)
ntERP_fmri = squeeze(mean(nontargetEEG_fmri, 1)); % ERP não-target

% ----- Deitado -----

tERP_lying = squeeze(mean(targetEEG_lying, 1)); % ERP target (média do sinal na dimensão trials)
ntERP_lying = squeeze(mean(nontargetEEG_lying, 1)); % ERP não-target

%% Correção da linha de base

% nobaseline = tERP_seated;
% 
% baseline = [-200;-100];
% 
% baseidx = dsearchn(time',baseline);
% 
% % ----- Sentado -----
% 
% for e=1:length(channel_labels)
%     tmn = mean(tERP_seated(baseidx(1):baseidx(2),e));   
%     ntmn = mean(ntERP_seated(baseidx(1):baseidx(2),e));   
% 
%     tERP_seated(:,e) = tERP_seated(:,e) - tmn;
%     ntERP_seated(:,e) = ntERP_seated(:,e) - ntmn;
% end
% 
% % ----- fMRI -----
% 
% for e=1:length(channel_labels)
%     tmn = mean(tERP_fmri(baseidx(1):baseidx(2),e));   
%     ntmn = mean(ntERP_fmri(baseidx(1):baseidx(2),e));   
% 
%     tERP_fmri(:,e) = tERP_fmri(:,e) - tmn;
%     ntERP_fmri(:,e) = ntERP_fmri(:,e) - ntmn;
% end
% 
% % ----- Deitado -----
% 
% for e=1:length(channel_labels)
%     tmn = mean(tERP_lying(baseidx(1):baseidx(2),e));   
%     ntmn = mean(ntERP_lying(baseidx(1):baseidx(2),e));   
% 
%     tERP_lying(:,e) = tERP_lying(:,e) - tmn;
%     ntERP_lying(:,e) = ntERP_lying(:,e) - ntmn;
% end

%% Baseline vs no baseline

% fig = fig + 1;
% figure(fig)
% hold on,
% xline(0, '--','LineWidth',1.5)
% plot(time,tERP_seated(:,8),'b','LineWidth',2);
% plot(time,nobaseline(:,8),'r','LineWidth',2);
% title("Baseline vs No baseline correction"), xlabel('Time (ms)'), ylabel('Amplitude (\muv)'), grid on
% ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[12 12 -8 -8],'black');
% set(ph0,'facealpha',.15,'edgecolor','none')
% legend('Stimuli Onset','Baseline corrected','No baseline corrected','Baseline period', 'AutoUpdate', 'off');
% yline(0, ':','LineWidth',1.5)
% hold off,

%% Single trial plot

fig = fig + 1;
figure(fig)
hold on, grid on,
xline(0, '--','LineWidth',1.5)
yline(0, ':','LineWidth',1.5)
title(nch_locs(8).labels)
plot(time,targetEEG_fmri(101,:,8),'b','LineWidth',2);

%% Gráficos - Sentado

% Definição do tempo dos picos e respetivos índices 

tN170 = [90; 130]; % ms
idxN170 = dsearchn(time',tN170);

tN2 = [180; 230];
idxN2 = dsearchn(time',tN2);

tP3 = [240; 340];
idxP3 = dsearchn(time',tP3);

%% Gráficos dos ERPs no tempo

i = 7;

save('tERP_seatedPO8.mat','tERP_seated')

clim = [floor(min(min(tERP_seated(:,i)))),ceil(max(max(tERP_seated(:,i))))];

aux = 0;

figure(fig)
hold on,
fig = fig + 1;
plot(time,tERP_seated(:,i),'b','LineWidth',2);
plot(time,ntERP_seated(:,i),'r','LineWidth',2);
xline(0, '--')    
ylim([(clim(1)+aux),(clim(2)-aux)])
xlim([min(time),max(time)])
title(channel_labels(i) + " - Seated"), xlabel('Time (ms)'), ylabel('Amplitude (\muv)'), grid on

if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
    c = 'blue';
    a = 'VPP';
    N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
else % N170
    c = 'cyan';
    a = 'N170';
    N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
end

N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];

% Estas patches precisam de ser ajustadas

ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
set(ph0,'facealpha',.1,'edgecolor','none')

ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
set(ph1,'facealpha',.1,'edgecolor','none')

ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
set(ph2,'facealpha',.1,'edgecolor','none')

ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
set(ph3,'facealpha',.1,'edgecolor','none')

legend('Target','Non-target','Stimulus','Baseline',a,'N2','AutoUpdate', 'off')

yline(0, '--')

hold off

%
clim = [floor(min(min(tERP_seatedsound(:,i)))),ceil(max(max(tERP_seatedsound(:,i))))];

aux = 0;

figure(fig)
hold on,
fig = fig + 1;
plot(time,tERP_seatedsound(:,i),'b','LineWidth',2);
plot(time,ntERP_seatedsound(:,i),'r','LineWidth',2);
xline(0, '--')    
ylim([(clim(1)+aux),(clim(2)-aux)])
xlim([min(time),max(time)])
title(channel_labels(i) + " - Seated w/ sound"), xlabel('Time (ms)'), ylabel('Amplitude (\muv)'), grid on

if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
    c = 'blue';
    a = 'VPP';
    N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
else % N170
    c = 'cyan';
    a = 'N170';
    N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
end

N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];

% Estas patches precisam de ser ajustadas

ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
set(ph0,'facealpha',.1,'edgecolor','none')

ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
set(ph1,'facealpha',.1,'edgecolor','none')

ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
set(ph2,'facealpha',.1,'edgecolor','none')

% ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
% set(ph3,'facealpha',.1,'edgecolor','none')

legend('Target','Non-target','Stimulus','Baseline',a,'N2','AutoUpdate', 'off')

yline(0, '--')

hold off

clim = [floor(min(min(tERP_lying(:,i)))),ceil(max(max(tERP_lying(:,i))))];

aux = 0;

figure(fig)
hold on,
fig = fig + 1;
plot(time,tERP_lying(:,i),'b','LineWidth',2);
plot(time,ntERP_lying(:,i),'r','LineWidth',2);
xline(0, '--')    
ylim([(clim(1)+aux),(clim(2)-aux)])
xlim([min(time),max(time)])
title(channel_labels(i) + " - Lying"), xlabel('Time (ms)'), ylabel('Amplitude (\muv)'), grid on

if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
    c = 'blue';
    a = 'VPP';
    N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
else % N170
    c = 'cyan';
    a = 'N170';
    N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
end

N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];

% Estas patches precisam de ser ajustadas

ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
set(ph0,'facealpha',.1,'edgecolor','none')

ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
set(ph1,'facealpha',.1,'edgecolor','none')

ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
set(ph2,'facealpha',.1,'edgecolor','none')

% ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
% set(ph3,'facealpha',.1,'edgecolor','none')

legend('Target','Non-target','Stimulus','Baseline',a,'N2','AutoUpdate', 'off')

yline(0, '--')

hold off

%
clim = [floor(min(min(tERP_fmri(:,i)))),ceil(max(max(tERP_fmri(:,i))))];

aux = 0;

figure(fig)
hold on,
fig = fig + 1;
plot(time,tERP_fmri(:,i),'b','LineWidth',2);
plot(time,ntERP_fmri(:,i),'r','LineWidth',2);
xline(0, '--')    
ylim([(clim(1)+aux),(clim(2)-aux)])
xlim([min(time),max(time)])
title(channel_labels(i) + " - Lying w/ fRMI"), xlabel('Time (ms)'), ylabel('Amplitude (\muv)'), grid on

if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
    c = 'blue';
    a = 'VPP';
    N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
else % N170
    c = 'cyan';
    a = 'N170';
    N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
end

N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];

% Estas patches precisam de ser ajustadas

ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
set(ph0,'facealpha',.1,'edgecolor','none')

ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
set(ph1,'facealpha',.1,'edgecolor','none')

ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
set(ph2,'facealpha',.1,'edgecolor','none')

% ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
% set(ph3,'facealpha',.1,'edgecolor','none')

legend('Target','Non-target','Stimulus','Baseline',a,'N2','AutoUpdate', 'off')

yline(0, '--')

hold off

% 
% % Mapas topográficos
% 
% figure(fig)
% fig = fig + 1;
% tlo1=tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
% title(tlo1, ["Sentado - Target (Em cima) vs Non-Target (Em baixo)";""], 'FontSize', 15, 'FontWeight','Bold');
% 
% ax1=nexttile(tlo1);
% topoplot(mean(tERP_seated(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax1,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)
% 
% ax2=nexttile(tlo1);
% topoplot(mean(tERP_seated(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax2,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)
% 
% ax3=nexttile(tlo1);
% topoplot(mean(tERP_seated(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax3,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)
% 
% ax4=nexttile(tlo1);
% topoplot(mean(ntERP_seated(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax4,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)
% 
% ax5=nexttile(tlo1);
% topoplot(mean(ntERP_seated(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax5,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)
% 
% ax6=nexttile(tlo1);
% topoplot(mean(ntERP_seated(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax6,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)
% 
% cbh2=colorbar(ax6);
% cbh2.Layout.Tile = 'east';
% cbh2.Label.String='\muV';
% chh2.Label.FontSize = 12;
% cbh2.Label.Rotation=0; 
% 
% %% Gráficos - fMRI
% 
% % Gráficos dos ERPs no tempo
% 
% clim = [floor(min(min(tERP_fmri))),ceil(max(max(tERP_fmri)))];
% 
% aux = 0;
% 
% for i=1:length(channel_labels)
%     figure(fig)
%     hold on,
%     fig = fig + 1;
%     plot(time,tERP_fmri(:,i),'b','LineWidth',2);
%     plot(time,ntERP_fmri(:,i),'r','LineWidth',2);
%     xline(0, '--')    
%     ylim([(clim(1)+aux),(clim(2)-aux)])
%     xlim([min(time),max(time)])
%     title(channel_labels(i) + " - EEG com fMRI"), xlabel('Tempo (ms)'), ylabel('Amplitude (\muv)'), grid on
% 
%     if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
%         c = 'blue';
%         a = 'VPP';
%         N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
%     else % N170
%         c = 'cyan';
%         a = 'N170';
%         N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
%     end
% 
%     N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
%     P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];
% 
%     % Estas patches precisam de ser ajustadas
% 
%     ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
%     set(ph0,'facealpha',.1,'edgecolor','none')
% 
%     ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
%     set(ph1,'facealpha',.1,'edgecolor','none')
% 
%     ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
%     set(ph2,'facealpha',.1,'edgecolor','none')
% 
%     ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
%     set(ph3,'facealpha',.1,'edgecolor','none')
% 
%     legend('Target','Non-target','Stimulus','Baseline',a, 'N2', 'P3','AutoUpdate', 'off')
% 
%     yline(0, '--')
%     hold off
% 
% end
% 
% % Mapas topográficos
% 
% figure(fig)
% fig = fig + 1;
% tlo1=tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
% title(tlo1, ["EEG com fMRI - Target (Em cima) vs Non-Target (Em baixo)";""], 'FontSize', 15, 'FontWeight','Bold');
% 
% ax1=nexttile(tlo1);
% topoplot(mean(tERP_fmri(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax1,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)
% 
% ax2=nexttile(tlo1);
% topoplot(mean(tERP_fmri(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax2,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)
% 
% ax3=nexttile(tlo1);
% topoplot(mean(tERP_fmri(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax3,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)
% 
% ax4=nexttile(tlo1);
% topoplot(mean(ntERP_fmri(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax4,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)
% 
% ax5=nexttile(tlo1);
% topoplot(mean(ntERP_fmri(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax5,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)
% 
% ax6=nexttile(tlo1);
% topoplot(mean(ntERP_fmri(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax6,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)
% 
% cbh2=colorbar(ax6);
% cbh2.Layout.Tile = 'east';
% cbh2.Label.String='\muV';
% chh2.Label.FontSize = 12;
% cbh2.Label.Rotation=0; 
% 
% %% Gráficos - Deitado
% 
% % Gráficos dos ERPs no tempo
% 
% clim = [floor(min(min(tERP_lying))),ceil(max(max(tERP_lying)))];
% 
% aux = 0;
% 
% for i=1:length(channel_labels)
%     figure(fig)
%     hold on,
%     fig = fig + 1;
%     plot(time,tERP_lying(:,i),'b','LineWidth',2);
%     plot(time,ntERP_lying(:,i),'r','LineWidth',2);
%     xline(0, '--')    
%     ylim([(clim(1)+aux),(clim(2)-aux)])
%     xlim([min(time),max(time)])
%     title(channel_labels(i) + " - Deitado"), xlabel('Tempo (ms)'), ylabel('Amplitude (\muv)'), grid on
% 
%     if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
%         c = 'blue';
%         a = 'VPP';
%         N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
%     else % N170
%         c = 'cyan';
%         a = 'N170';
%         N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
%     end
% 
%     N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
%     P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];
% 
%     % Estas patches precisam de ser ajustadas
% 
%     ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
%     set(ph0,'facealpha',.1,'edgecolor','none')
% 
%     ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
%     set(ph1,'facealpha',.1,'edgecolor','none')
% 
%     ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
%     set(ph2,'facealpha',.1,'edgecolor','none')
% 
%     ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
%     set(ph3,'facealpha',.1,'edgecolor','none')
% 
%     legend('Target','Non-target','Stimulus','Baseline',a, 'N2', 'P3','AutoUpdate', 'off')
% 
%     yline(0, '--')
%     hold off
% 
% end
% 
% % Mapas topográficos
% 
% figure(fig)
% fig = fig + 1;
% tlo1=tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
% title(tlo1, ["Deitado - Target (Em cima) vs Non-Target (Em baixo)";""], 'FontSize', 15, 'FontWeight','Bold');
% 
% ax1=nexttile(tlo1);
% topoplot(mean(tERP_lying(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax1,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)
% 
% ax2=nexttile(tlo1);
% topoplot(mean(tERP_lying(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax2,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)
% 
% ax3=nexttile(tlo1);
% topoplot(mean(tERP_lying(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax3,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)
% 
% ax4=nexttile(tlo1);
% topoplot(mean(ntERP_lying(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax4,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)
% 
% ax5=nexttile(tlo1);
% topoplot(mean(ntERP_lying(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax5,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)
% 
% ax6=nexttile(tlo1);
% topoplot(mean(ntERP_lying(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
% title(ax6,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)
% 
% cbh2=colorbar(ax6);
% cbh2.Layout.Tile = 'east';
% cbh2.Label.String='\muV';
% chh2.Label.FontSize = 12;
% cbh2.Label.Rotation=0; 

%% Comparação entre condições - Target

% Gráficos dos ERPs no tempo

clim = [floor(min(min(tERP_seated))),ceil(max(max(tERP_seated))),
        floor(min(min(tERP_fmri))),ceil(max(max(tERP_fmri))),
        floor(min(min(tERP_lying))),ceil(max(max(tERP_lying)))];

clim = [min(min(clim)),max(max(clim))];

aux = 0;

for i=1:length(channel_labels)
    figure(fig)
    hold on,
    fig = fig + 1;
    plot(time,tERP_seated(:,i),'b','LineWidth',2);
    plot(time,tERP_seatedsound(:,i),'k','LineWidth',2);
    plot(time,tERP_lying(:,i),'r','LineWidth',2);
    plot(time,tERP_fmri(:,i),'g','LineWidth',2);
    xline(0, '--','LineWidth',1.5)    
    ylim([(clim(1)+aux),(clim(2)-aux)])
    xlim([min(time),max(time)])
    title(channel_labels(i) + " - Target"), xlabel('Time (ms)'), ylabel('Amplitude (\muv)'), grid on
    
    if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
        c = 'blue';
        a = 'VPP';
        N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
    else % N170
        c = 'cyan';
        a = 'N170';
        N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
    end
   
    N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
    P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];

    % Estas patches precisam de ser ajustadas

    ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
    set(ph0,'facealpha',.1,'edgecolor','none')

    % ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
    % set(ph1,'facealpha',.1,'edgecolor','none')
    % 
    % ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
    % set(ph2,'facealpha',.1,'edgecolor','none')
    
    % ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
    % set(ph3,'facealpha',.1,'edgecolor','none')

    legend('Seated','Seated w/ fMRI', 'Lying', 'Lying w/ fMRI','Stimulus','Baseline', 'P3','AutoUpdate', 'off')
    
    yline(0, ':','LineWidth',1.5)
    hold off
    
end

%% Mapas topográficos - Target

figure(fig)
fig = fig + 1;
tlo1=tiledlayout(4,3,'TileSpacing','none','Padding','none');
title(tlo1, ["Target";""], 'FontSize', 15, 'FontWeight','Bold');

n_cont = 10;

ax1=nexttile(tlo1);
topoplot(mean(tERP_seated(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');
title(ax1,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)

ax2=nexttile(tlo1);
topoplot(mean(tERP_seated(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');
title(ax2,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)

ax3=nexttile(tlo1);
topoplot(mean(tERP_seated(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');
title(ax3,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)

ax14=nexttile(tlo1);
topoplot(mean(tERP_seatedsound(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax15=nexttile(tlo1);
topoplot(mean(tERP_seatedsound(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax16=nexttile(tlo1);
topoplot(mean(tERP_seatedsound(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax4=nexttile(tlo1);
topoplot(mean(tERP_lying(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax4=nexttile(tlo1);
topoplot(mean(tERP_lying(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax6=nexttile(tlo1);
topoplot(mean(tERP_lying(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax7=nexttile(tlo1);
topoplot(mean(tERP_fmri(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax8=nexttile(tlo1);
topoplot(mean(tERP_fmri(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax9=nexttile(tlo1);
topoplot(mean(tERP_fmri(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

cbh2=colorbar(ax9);
cbh2.Layout.Tile = 'east';
cbh2.Label.String='\muV';
chh2.Label.FontSize = 16;
cbh2.Label.Rotation=0;

%% Comparação entre condições - Non-Target

% Gráficos dos ERPs no tempo

aux = 0;

for i=1:length(channel_labels)
    figure(fig)
    hold on,
    fig = fig + 1;
    plot(time,ntERP_seated(:,i),'b','LineWidth',2);
    plot(time,ntERP_lying(:,i),'r','LineWidth',2);
    plot(time,ntERP_fmri(:,i),'g','LineWidth',2);
    xline(0, '--')    
    ylim([(clim(1)+aux),(clim(2)-aux)])
    xlim([min(time),max(time)])
    title(channel_labels(i) + " - Non-Target"), xlabel('Tempo (ms)'), ylabel('Amplitude (\muv)'), grid on
    
    if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
        c = 'blue';
        a = 'VPP';
        N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
    else % N170
        c = 'cyan';
        a = 'N170';
        N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
    end
   
    N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
    P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];

    % Estas patches precisam de ser ajustadas

    ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
    set(ph0,'facealpha',.1,'edgecolor','none')

    ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
    set(ph1,'facealpha',.1,'edgecolor','none')

    ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
    set(ph2,'facealpha',.1,'edgecolor','none')
    
    ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
    set(ph3,'facealpha',.1,'edgecolor','none')

    legend('EEG Seated','EEG Lying','EEG w/ fMRI','Stimulus','Baseline',a, 'N2', 'P3','AutoUpdate', 'off')
    
    yline(0, ':')
    hold off
    
end

%% Mapas topográficos - Non-target

clim = [floor(min(min(tERP_seated))),ceil(max(max(tERP_seated))),
        floor(min(min(tERP_fmri))),ceil(max(max(tERP_fmri))),
        floor(min(min(tERP_lying))),ceil(max(max(tERP_lying)))];

clim = [min(min(clim)),max(max(clim))];

figure(fig)
fig = fig + 1;
tlo1=tiledlayout(4,3,'TileSpacing','none','Padding','none');
title(tlo1, ["NonTarget";""], 'FontSize', 15, 'FontWeight','Bold');

ax1=nexttile(tlo1);
topoplot(mean(ntERP_seated(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');
title(ax1,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)

ax2=nexttile(tlo1);
topoplot(mean(ntERP_seated(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');
title(ax2,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)

ax3=nexttile(tlo1);
topoplot(mean(ntERP_seated(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');
title(ax3,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)

ax14=nexttile(tlo1);
topoplot(mean(ntERP_seatedsound(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax15=nexttile(tlo1);
topoplot(mean(ntERP_seatedsound(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax16=nexttile(tlo1);
topoplot(mean(ntERP_seatedsound(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax4=nexttile(tlo1);
topoplot(mean(ntERP_lying(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax5=nexttile(tlo1);
topoplot(mean(ntERP_lying(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax6=nexttile(tlo1);
topoplot(mean(ntERP_lying(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax7=nexttile(tlo1);
topoplot(mean(ntERP_fmri(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax8=nexttile(tlo1);
topoplot(mean(ntERP_fmri(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

ax9=nexttile(tlo1);
topoplot(mean(ntERP_fmri(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',n_cont,'electrodes','on');

cbh2=colorbar(ax9);
cbh2.Layout.Tile = 'east';
cbh2.Label.String='\muV';
chh2.Label.FontSize = 16;
cbh2.Label.Rotation=0;

%% Comparação entre condições - Target vs Non-Target

% Gráficos dos ERPs no tempo

% aux = 0;
% 
% for i=1:length(channel_labels)
%     figure(fig)
%     hold on,
%     fig = fig + 1;
%     plot(time,tERP_seated(:,i),'b','LineWidth',2);
%     plot(time,tERP_lying(:,i),'r','LineWidth',2);
%     plot(time,tERP_fmri(:,i),'g','LineWidth',2);
%     plot(time,ntERP_seated(:,i),'b--','LineWidth',2);
%     plot(time,ntERP_lying(:,i),'r--','LineWidth',2);
%     plot(time,ntERP_fmri(:,i),'g--','LineWidth',2);
%     xline(0, '--')    
%     ylim([(clim(1)+aux),(clim(2)-aux)])
%     xlim([min(time),max(time)])
%     title(channel_labels(i) + " - Target vs Non-Target"), xlabel('Tempo (ms)'), ylabel('Amplitude (\muv)'), grid on
% 
%     if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
%         c = 'blue';
%         a = 'VPP';
%         N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
%     else % N170
%         c = 'cyan';
%         a = 'N170';
%         N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
%     end
% 
%     N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
%     P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];
% 
%     % Estas patches precisam de ser ajustadas
% 
%     ph0 = patch([baseline(1) baseline(2) baseline(2) baseline(1)],[clim(1) clim(1) clim(2) clim(2)],'black');
%     set(ph0,'facealpha',.1,'edgecolor','none')
% 
%     ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
%     set(ph1,'facealpha',.1,'edgecolor','none')
% 
%     ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
%     set(ph2,'facealpha',.1,'edgecolor','none')
% 
%     ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
%     set(ph3,'facealpha',.1,'edgecolor','none')
% 
%     legend('EEG Sentado - Target','EEG Deitado - Target','EEG c/ fMRI - Target','EEG Sentado - NonTarget','EEG Deitado - NonTarget','EEG c/ fMRI - NonTarget', ...      
%            'Stimulus','Baseline',a, 'N2', 'P3','AutoUpdate', 'off')
% 
%     yline(0, ':')
%     hold off
% 
% end

% %% 2D movie of topographic ERPs (Target - Seated)
% 
% figure('units','pixels','position',[0 0 1920 1080]); 
% set(gcf, 'WindowState', 'maximized');
% 
% targ_seated = eegmovie(tERP_seated', fs, nch_locs ,...
%      'time','on','startsec',-0.200 , 'topoplotopt', ...
%      {'numcontour' n_cont},'minmax',clim,'title','EEG seated');
% 
% seemovie(targ_seated,0,clim);
% 
% %% 2D movie of topographic ERPs (Target - Lying)
% 
% figure('units','pixels','position',[0 0 1920 1080]); 
% set(gcf, 'WindowState', 'maximized');
% 
% targ_lying = eegmovie(tERP_lying', fs, nch_locs ,...
%      'time','on','startsec',-0.200 , 'topoplotopt', ...
%      {'numcontour' n_cont},'minmax',clim,'title','EEG lying');
% 
% seemovie(targ_lying,0,clim);
% 
% %% 2D movie of topographic ERPs (Target - fMRI)
% 
% figure('units','pixels','position',[0 0 1920 1080]); 
% set(gcf, 'WindowState', 'maximized');
% 
% targ_fmri = eegmovie(tERP_fmri', fs, nch_locs ,...
%      'time','on','startsec',-0.200 , 'topoplotopt', ...
%      {'numcontour' n_cont},'minmax',clim,'title','EEG w/ fMRI');
% 
% seemovie(targ_fmri,0,clim);
% 
% %% 2D movie of topographic ERPs (NonTarget - Seated)
% 
% figure('units','pixels','position',[0 0 1920 1080]); 
% set(gcf, 'WindowState', 'maximized');
% 
% nontarg_seated = eegmovie(ntERP_seated', fs, nch_locs ,...
%      'time','on','startsec',-0.200 , 'topoplotopt', ...
%      {'numcontour' n_cont},'minmax',clim,'title','EEG seated');
% 
% seemovie(nontarg_seated,0,clim);
% 
% %% 2D movie of topographic ERPs (NonTarget - Lying)
% 
% figure('units','pixels','position',[0 0 1920 1080]); 
% set(gcf, 'WindowState', 'maximized');
% 
% nontarg_lying = eegmovie(ntERP_lying', fs, nch_locs ,...
%      'time','on','startsec',-0.200 , 'topoplotopt', ...
%      {'numcontour' n_cont},'minmax',clim,'title','EEG lying');
% 
% seemovie(nontarg_lying,0,clim);
% 
% %% 2D movie of topographic ERPs (NonTarget - fMRI)
% 
% figure('units','pixels','position',[0 0 1920 1080]); 
% set(gcf, 'WindowState', 'maximized');
% 
% nontarg_fmri = eegmovie(ntERP_fmri', fs, nch_locs ,...
%      'time','on','startsec',-0.200 , 'topoplotopt', ...
%      {'numcontour' n_cont},'minmax',clim,'title','EEG w/ fMRI');
% 
% seemovie(nontarg_fmri,0,clim);