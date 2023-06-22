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
signal_seated = [];

for i=1:length(seg_target_seated)
    signal_seated = [signal_seated; filt_data_seated((seg_target_seated(i)+idx_bf):(seg_target_seated(i)+idx_af),7)'];
end

% -------------------
% Target - fMRI

%[targetEEG_fmri, nontargetEEG_fmri] = deal(zeros(length(seg_target_fmri), length(time), length(channel_labels))); % Pré-alocação dos dados
signal_fMRI = [];

for i=1:length(seg_target_fmri)
    signal_fMRI = [signal_fMRI; filt_data_fmri((seg_target_fmri(i)+idx_bf):(seg_target_fmri(i)+idx_af),7)'];
end

% -------------------
% Target - Deitado

%[targetEEG_lying, nontargetEEG_lying] = deal(zeros(length(seg_target_lying), length(time), length(channel_labels))); % Pré-alocação dos dados
signal_lying = [];

for i=1:length(seg_target_lying)
    signal_lying = [signal_lying; filt_data_lying((seg_target_lying(i)+idx_bf):(seg_target_lying(i)+idx_af),7)'];
end

%% Correção da linha de base

baseline = [-200;0];

baseidx = dsearchn(time',baseline);

% ----- Sentado -----

for e=1:length(signal_seated)
    tmn = mean(signal_seated(e,baseidx(1):baseidx(2)));   
    signal_seated(e,:) = signal_seated(e,:) - tmn;
end

% ----- fMRI -----

for e=1:length(signal_fMRI)
    tmn = mean(signal_fMRI(e,baseidx(1):baseidx(2)));   
    signal_fMRI(e,:) = signal_fMRI(e,:) - tmn;
end

% ----- Deitado -----

for e=1:length(signal_lying)
    tmn = mean(signal_lying(e,baseidx(1):baseidx(2)));   
    signal_lying(e,:) = signal_lying(e,:) - tmn;
end


%% ANOVA P300

startTime = 0.240 - bf; startIndex = round(startTime * fs);
endTime = 0.340 - bf; endIndex = round(endTime * fs);

p300Interval_seated = signal_seated(:,startIndex:endIndex);
p300Interval_fMRI = signal_fMRI(:,startIndex:endIndex); 
p300Interval_lying = signal_lying(:,startIndex:endIndex);

[~, maxIndex_seated_p300] = max(p300Interval_seated');
[~, maxIndex_fMRI_p300] = max(p300Interval_fMRI');
[~, maxIndex_lying_p300] = max(p300Interval_lying');

for i = 1:length(p300Interval_seated)
    amplitude_seated_p300(i) = p300Interval_seated(i, maxIndex_seated_p300(i));
end

for i = 1:length(p300Interval_fMRI)
    amplitude_fMRI_p300(i) = p300Interval_fMRI(i, maxIndex_fMRI_p300(i));
end

for i = 1:length(p300Interval_lying)
    amplitude_lying_p300(i) = p300Interval_lying(i, maxIndex_lying_p300(i));
end

num_p300 = min([length(amplitude_seated_p300), length(amplitude_fMRI_p300), length(amplitude_lying_p300)]);

amplitude_seated_p300 = amplitude_seated_p300(1:num_p300);
amplitude_fMRI_p300 = amplitude_fMRI_p300(1:num_p300);
amplitude_lying_p300 = amplitude_lying_p300(1:num_p300);

anova_amp_p300 = [amplitude_seated_p300', amplitude_fMRI_p300', amplitude_lying_p300'];
[pValue_amp_p300, ~, stats_amp_p300] = anova1(anova_amp_p300);
fprintf('amplitude P300 P-value: %.4f\n', pValue_amp_p300);

latency_seated_p300 = (maxIndex_seated_p300./fs)+bf;
latency_fMRI_p300 = (maxIndex_fMRI_p300./fs)+bf;
latency_lying_p300 = (maxIndex_lying_p300./fs)+bf;

latency_seated_p300 = latency_seated_p300(1:num_p300);
latency_fMRI_p300 = latency_fMRI_p300(1:num_p300);
latency_lying_p300 = latency_lying_p300(1:num_p300);

anova_lat_p300 = [latency_seated_p300', latency_fMRI_p300', latency_lying_p300'];
[pValue_lat_p300, ~, stats_lat_p300] = anova1(anova_lat_p300);
fprintf('latency P300 P-value: %.4f\n', pValue_lat_p300);

%% ANOVA N170

startTime = 0.090 - bf; startIndex = round(startTime * fs);
endTime = 0.120 - bf; endIndex = round(endTime * fs);

n170Interval_seated = signal_seated(:,startIndex:endIndex);
n170Interval_fMRI = signal_fMRI(:,startIndex:endIndex); 
n170Interval_lying = signal_lying(:,startIndex:endIndex);

[~, minIndex_seated_n170] = min(n170Interval_seated');
[~, minIndex_fMRI_n170] = min(n170Interval_fMRI');
[~, minIndex_lying_n170] = min(n170Interval_lying');

for i = 1:length(n170Interval_seated)
    amplitude_seated_n170(i) = n170Interval_seated(i, minIndex_seated_n170(i));
end

for i = 1:length(n170Interval_fMRI)
    amplitude_fMRI_n170(i) = n170Interval_fMRI(i, minIndex_fMRI_n170(i));
end

for i = 1:length(n170Interval_lying)
    amplitude_lying_n170(i) = n170Interval_lying(i, minIndex_lying_n170(i));
end

num_n170 = min([length(amplitude_seated_n170), length(amplitude_fMRI_n170), length(amplitude_lying_n170)]);

amplitude_seated_n170 = amplitude_seated_n170(1:num_n170);
amplitude_fMRI_n170= amplitude_fMRI_n170(1:num_n170);
amplitude_lying_n170 = amplitude_lying_n170(1:num_n170);

anova_amp_n170 = [amplitude_seated_n170', amplitude_fMRI_n170', amplitude_lying_n170'];
[pValue_amp_n170, ~, stats_amp_n170] = anova1(anova_amp_n170);
fprintf('amplitude n170 P-value: %.4f\n', pValue_amp_n170);

latency_seated_n170 = (minIndex_seated_n170./fs)+bf;
latency_fMRI_n170 = (minIndex_fMRI_n170./fs)+bf;
latency_lying_n170 = (minIndex_lying_n170./fs)+bf;

latency_seated_n170 = latency_seated_n170(1:num_n170);
latency_fMRI_n170 = latency_fMRI_n170(1:num_n170);
latency_lying_n170 = latency_lying_n170(1:num_n170);

anova_lat_n170 = [latency_seated_n170', latency_fMRI_n170', latency_lying_n170'];
[pValue_amp_n170, ~, stats_lat_n170] = anova1(anova_lat_n170);
fprintf('latency n170 P-value: %.4f\n', pValue_amp_n170);
% fprintf('F-statistic: %.4f\n', stats{2, 5});
% fprintf('Degrees of freedom (between groups): %d\n', stats{2, 3});
% fprintf('Degrees of freedom (within groups): %d\n', stats{3, 3});

%% ANOVA N200

startTime = 0.180 - bf; startIndex = round(startTime * fs);
endTime = 0.220 - bf; endIndex = round(endTime * fs);

n200Interval_seated = signal_seated(:,startIndex:endIndex);
n200Interval_fMRI = signal_fMRI(:,startIndex:endIndex); 
n200Interval_lying = signal_lying(:,startIndex:endIndex);

[~, minIndex_seated_n200] = min(n200Interval_seated');
[~, minIndex_fMRI_n200] = min(n200Interval_fMRI');
[~, minIndex_lying_n200] = min(n200Interval_lying');

for i = 1:length(n200Interval_seated)
    amplitude_seated_n200(i) = n200Interval_seated(i, minIndex_seated_n200(i));
end

for i = 1:length(n200Interval_fMRI)
    amplitude_fMRI_n200(i) = n200Interval_fMRI(i, minIndex_fMRI_n200(i));
end

for i = 1:length(n200Interval_lying)
    amplitude_lying_n200(i) = n200Interval_lying(i, minIndex_lying_n200(i));
end

num_n200 = min([length(amplitude_seated_n200), length(amplitude_fMRI_n200), length(amplitude_lying_n200)]);

amplitude_seated_n200 = amplitude_seated_n200(1:num_n200);
amplitude_fMRI_n200= amplitude_fMRI_n200(1:num_n200);
amplitude_lying_n200 = amplitude_lying_n200(1:num_n200);

anova_amp_n200 = [amplitude_seated_n200', amplitude_fMRI_n200', amplitude_lying_n200'];
[pValue_amp_n200, ~, stats_amp_n200] = anova1(anova_amp_n200);
fprintf('amplitude n200 P-value: %.4f\n', pValue_amp_n200);

latency_seated_n200 = (minIndex_seated_n200./fs)+bf;
latency_fMRI_n200 = (minIndex_fMRI_n200./fs)+bf;
latency_lying_n200 = (minIndex_lying_n200./fs)+bf;

latency_seated_n200 = latency_seated_n200(1:num_n170);
latency_fMRI_n200 = latency_fMRI_n200(1:num_n170);
latency_lying_n200 = latency_lying_n200(1:num_n170);

anova_lat_n200 = [latency_seated_n200', latency_fMRI_n200', latency_lying_n200'];
[pValue_amp_n200, ~, stats_lat_n200] = anova1(anova_lat_n200);
fprintf('latency n200 P-value: %.4f\n', pValue_amp_n200);
% fprintf('F-statistic: %.4f\n', stats{2, 5});
% fprintf('Degrees of freedom (between groups): %d\n', stats{2, 3});
% fprintf('Degrees of freedom (within groups): %d\n', stats{3, 3});