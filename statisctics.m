%% Importar os dados e inicializar variáveis

close all;
clc, clear;

% EEGLAB functions documentation: https://sccn.ucsd.edu/~arno/eeglab/auto/indexfunc.html

addpath(genpath('eeglab2023.0')); % Path do EEGLAB 

data_seated = readtable("S3_25_Male.csv");
data_seated = table2array(data_seated);

data_fmri = readtable("S2_25_Male_2023-06-02_16-21.csv");
data_fmri = table2array(data_fmri);

data_lying = readtable("S22_25_Male_2023-06-02_17-03.csv");
data_lying = table2array(data_lying);

fs = 125;
channel_labels = ["C3";"CP3";"P3";"P7";"PO7";"PO3";"Pz";"CPz";"Fz";"Cz";"C4";"CP4";"P4";"P8";"PO8";"PO4"];

%% Dados da fotoresistência 

photoVector_seated=(data_seated(:,end));

t_seated = ((0:length(photoVector_seated)-1)*1/fs)/60;
n_pnts = 20;
filt_photoVector_seated = movmedian(photoVector_seated,20);

photoVector_fmri=(data_fmri(:,end));

t_fmri = ((0:length(photoVector_fmri)-1)*1/fs)/60;
filt_photoVector_fmri = movmedian(photoVector_fmri,20);

photoVector_lying=(data_lying(:,end));

t_lying = ((0:length(photoVector_lying)-1)*1/fs)/60;
filt_photoVector_lying = movmedian(photoVector_lying,20);

%% Re-referenciação do EEG 

eeg_data_seated = data_seated(:,1:end-1);
eeg_data_fmri = data_fmri(:,1:end-1);
eeg_data_lying = data_lying(:,1:end-1);

%% Pré-processamento do EEG

% Band-pass
bpFilt = designfilt('bandpassfir','FilterOrder',3000, ...
         'CutoffFrequency1',0.1,'CutoffFrequency2',30, ...
         'SampleRate',fs); 

filt_data_seated = filtfilt(bpFilt,eeg_data_seated);
filt_data_fmri = filtfilt(bpFilt,eeg_data_fmri);
filt_data_lying = filtfilt(bpFilt,eeg_data_lying);

%% Retirar sub-harmónica dos 25 Hz

% Band-stop
bsFilt = designfilt('bandstopfir','FilterOrder',1000, ...
         'CutoffFrequency1',24.7,'CutoffFrequency2',25.3, ...
         'SampleRate',fs);
%fvtool(bsFilt)

filt_data_seated = filtfilt(bsFilt,filt_data_seated);
filt_data_fmri = filtfilt(bsFilt,filt_data_fmri);
filt_data_lying = filtfilt(bsFilt,filt_data_lying);

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

    elseif photoresistor >= 230 % Possible rising edge
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
    if(min(filt_photoVector_seated(falling_edges(i):rising_edges(i))) < 125)
        seg_nontarget_seated = [seg_nontarget_seated, falling_edges(i)];
    else
        seg_target_seated = [seg_target_seated, falling_edges(i)];
    end
end

seg_target_seated = unique(seg_target_seated);
seg_nontarget_seated = unique(seg_nontarget_seated);

% ----- fMRI -----

last_low = true; 
last_high = false;

falling_edges = [];
rising_edges = [];

for i = 1:length(filt_photoVector_fmri)
    photoresistor = filt_photoVector_fmri(i);
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

rising_edges(1:11) = []; 
falling_edges(1:10) = []; 
falling_edges(end) = [];

seg_target_fmri = [];
seg_nontarget_fmri = [];

for i = 1:length(falling_edges) % Verificar quais blocos pertecem aos targets e não targets
    if(min(filt_photoVector_fmri(falling_edges(i):rising_edges(i))) < 125)
        seg_nontarget_fmri = [seg_nontarget_fmri, falling_edges(i)];
    else
        seg_target_fmri = [seg_target_fmri, falling_edges(i)];
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
    if(min(filt_photoVector_lying(falling_edges(i):rising_edges(i))) < 125)
        seg_nontarget_lying = [seg_nontarget_lying, falling_edges(i)];
    else
        seg_target_lying = [seg_target_lying, falling_edges(i)];
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

%targetEEG_seated = zeros(length(seg_target_seated), length(time), length(channel_labels)); % Pré-alocação dos dados
p300_seated = [];

for i=1:length(seg_target_seated)
    p300_seated = [p300_seated; filt_data_seated((seg_target_seated(i)+idx_bf):(seg_target_seated(i)+idx_af),7)'];
end

% -------------------
% Target - fMRI

%[targetEEG_fmri, nontargetEEG_fmri] = deal(zeros(length(seg_target_fmri), length(time), length(channel_labels))); % Pré-alocação dos dados
p300_fMRI = [];

for i=1:length(seg_target_fmri)
    p300_fMRI = [p300_fMRI; filt_data_fmri((seg_target_fmri(i)+idx_bf):(seg_target_fmri(i)+idx_af),7)'];
end

% -------------------
% Target - Deitado

%[targetEEG_lying, nontargetEEG_lying] = deal(zeros(length(seg_target_lying), length(time), length(channel_labels))); % Pré-alocação dos dados
p300_lying = [];

for i=1:length(seg_target_lying)
    p300_lying = [p300_lying; filt_data_lying((seg_target_lying(i)+idx_bf):(seg_target_lying(i)+idx_af),7)'];
end

%% Correção da linha de base

baseline = [-200;0];

baseidx = dsearchn(time',baseline);

% ----- Sentado -----

for e=1:length(p300_seated)
    tmn = mean(p300_seated(e,baseidx(1):baseidx(2)));   
    p300_seated(e,:) = p300_seated(e,:) - tmn;
end

% ----- fMRI -----

for e=1:length(p300_fMRI)
    tmn = mean(p300_fMRI(e,baseidx(1):baseidx(2)));   
    p300_fMRI(e,:) = p300_fMRI(e,:) - tmn;
end

% ----- Deitado -----

for e=1:length(p300_lying)
    tmn = mean(p300_lying(e,baseidx(1):baseidx(2)));   
    p300_lying(e,:) = p300_lying(e,:) - tmn;
end


%% ANOVA P300 amp

startTime = 0.250 - bf; startIndex = round(startTime * fs);
endTime = 0.500 - bf; endIndex = round(endTime * fs);

p300Interval_seated = p300_seated(:,startIndex:endIndex);
p300Interval_fMRI = p300_fMRI(:,startIndex:endIndex); 
p300Interval_lying = p300_lying(:,startIndex:endIndex);

[~, maxIndex_seated] = max(p300Interval_seated');
[~, maxIndex_fMRI] = max(p300Interval_fMRI');
[~, maxIndex_lying] = max(p300Interval_lying');

for i = 1:length(p300Interval_seated)
    amplitude_seated(i) = p300Interval_seated(i, maxIndex_seated(i));
end

for i = 1:length(p300Interval_fMRI)
    amplitude_fMRI(i) = p300Interval_fMRI(i, maxIndex_fMRI(i));
end

for i = 1:length(p300Interval_lying)
    amplitude_lying(i) = p300Interval_lying(i, maxIndex_lying(i));
end

num = min([length(amplitude_seated), length(amplitude_fMRI), length(amplitude_lying)]);

amplitude_seated = amplitude_seated(1:num);
amplitude_fMRI = amplitude_fMRI(1:num);
amplitude_lying = amplitude_lying(1:num);

anova_amp = [amplitude_seated', amplitude_fMRI', amplitude_lying'];
[pValue, ~, stats_amp] = anova1(anova_amp);
fprintf('amplitude P-value: %.4f\n', pValue);
% fprintf('F-statistic: %.4f\n', stats{2, 5});
% fprintf('Degrees of freedom (between groups): %d\n', stats{2, 3});
% fprintf('Degrees of freedom (within groups): %d\n', stats{3, 3});

%% Get mean and sd P300 latency

latency_seated = (maxIndex_seated./fs)+bf;
latency_fMRI = (maxIndex_fMRI./fs)+bf;
latency_lying = (maxIndex_lying./fs)+bf;

latency_seated = latency_seated(1:num);
latency_fMRI = latency_fMRI(1:num);
latency_lying = latency_lying(1:num);

anova_lat = [latency_seated', latency_fMRI', latency_lying'];
[pValue, ~, stats_lat] = anova1(anova_lat);
fprintf('latency P-value: %.4f\n', pValue);
% fprintf('F-statistic: %.4f\n', stats{2, 5});
% fprintf('Degrees of freedom (between groups): %d\n', stats{2, 3});
% fprintf('Degrees of freedom (within groups): %d\n', stats{3, 3});