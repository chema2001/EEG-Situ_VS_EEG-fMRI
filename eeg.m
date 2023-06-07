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

addpath(genpath('eeglab2023.0')); % Path do EEGLAB
ch_locs=readlocs('BioSemi64.loc'); % Ficheiro com localização dos elétrodos

data_seated = readtable("S3_25_Male.csv");
data_seated = table2array(data_seated);

data_fmri = readtable("S2_25_Male_2023-06-02_16-21.csv");
data_fmri = table2array(data_fmri);

data_lying = readtable("S22_25_Male_2023-06-02_17-03.csv");
data_lying = table2array(data_lying);

fs = 125;
channel_labels = ["C3";"CP3";"P3";"PO3";"P7";"PO7";"Fz";"Cz";"CPz";"Pz";"C4";"CP4";"P4";"PO4";"P8";"PO8"];

idx = zeros(length(channel_labels),1);

for i = 1:length(ch_locs)
    cur_lb = string(ch_locs(i).labels);
    if (nonzeros(channel_labels == cur_lb))
        idx(i) = 1;
    end
end

nch_locs = ch_locs(logical(idx)); % Elétrodos usados


%% Verificar o efeito dos filtros na fase dos sinais

% t = (-3:1/fs:3)';
% s = 5*square(2*pi*t); % Sinal com a mesma frequência que o sinal da fotoresistência
% 
% sqwave = awgn(s, 20, 'measured'); % Gaussian white noise
% 
% plot(t, [sqwave, movmean(sqwave,10), movmedian(sqwave,20)]),
% title('Verificar o efeito dos filtros na fase dos sinais'), 
% xlabel('Tempo (s)'),
% legend('Original','Média flutuante','Mediana flutuante');

% Resumindo, a mediana flutuante não distorce tanto as fases do sinal,
% principalmente no início das rising edges e falling edges. Isto é
% extremamente importante em análises como as dos ERPs pois é preciso
% saber com exatidão o instante em que o estímulo é apresentado.

%% Dados da fotoresistência 

photoVector_seated=(data_seated(:,end));

t_seated = ((0:length(photoVector_seated)-1)*1/fs)/60;
n_pnts = 20;
filt_photoVector_seated = movmedian(photoVector_seated,20);

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

fig = 1;

%figure (fig); plot(t,photoVector), title('Fotoresistência'), xlabel('Tempo (min)'); % Sem filtros
%fig = fig + 1;

% figure (fig); plot(t_seated,filt_photoVector_seated), title("Fotoresistência (Mediana flutuante com "+n_pnts+" pontos) - EEG Sentado"), xlabel('Tempo (min)');
% fig = fig + 1;

% figure (fig); plot(t,lpFilt_phVec), title("Fotoresistência (Filtro passa-baixo)"), xlabel('Tempo (min)');
% fig = fig + 1;

% figure (fig); plot(t_fmri,filt_photoVector_fmri), title("Fotoresistência (Mediana flutuante com "+n_pnts+" pontos) - EEG com fMRI"), xlabel('Tempo (min)');
% fig = fig + 1;

% figure (fig); plot(t_lying,filt_photoVector_lying), title("Fotoresistência (Mediana flutuante com "+n_pnts+" pontos) - EEG Deitado"), xlabel('Tempo (min)');
% fig = fig + 1;

%% Re-referenciação do EEG 

eeg_data_seated = data_seated(:,1:end-1);
eeg_data_fmri = data_fmri(:,1:end-1);
eeg_data_lying = data_lying(:,1:end-1);

%car = bsxfun(@minus, eeg_data, mean(eeg_data,2));
%lap = laplacian_perrinX(eeg_data',[nch_locs.X],[nch_locs.Y],[nch_locs.Z]);

%% Pré-processamento do EEG

% Filtragem
% De acordo com o artigo com a literatura dos ERPs, um
% filtro passa-banda entre 0.5 e 30/40 Hz será usado.

% Band-pass
bpFilt = designfilt('bandpassfir','FilterOrder',1000, ...
         'CutoffFrequency1',0.5,'CutoffFrequency2',30, ...
         'SampleRate',fs); 
%fvtool(bpFilt) % Visualizar resposta em frequência e fase do filtro

%filt_data = filtfilt(bpFilt,car);
%filt_data = filtfilt(bpFilt,lap');

filt_data_seated = filtfilt(bpFilt,eeg_data_seated);
filt_data_fmri = filtfilt(bpFilt,eeg_data_fmri);
filt_data_lying = filtfilt(bpFilt,eeg_data_lying);

%% FFT do sinal 

% Verificar se o sinal é um EEG. Para isso, o sinal no domíno das
% frequências tem que seguir a forma 1/f.
% Isto após retirar as componentes DC e dos 50Hz da rede elétrica.

% ----- Sentado -----

T = 1/fs; % Sampling period
L = length(filt_data_seated(:,8)); % Length of signal
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

bf = -0.100; % 100 ms antes do estímulo
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

%% Computação dos ERPs

% ----- Sentado -----

tERP_seated = squeeze(mean(targetEEG_seated, 1)); % ERP target (média do sinal na dimensão trials)
ntERP_seated = squeeze(mean(nontargetEEG_seated, 1)); % ERP não-target

% ----- fMRI -----

tERP_fmri = squeeze(mean(targetEEG_fmri, 1)); % ERP target (média do sinal na dimensão trials)
ntERP_fmri = squeeze(mean(nontargetEEG_fmri, 1)); % ERP não-target

% ----- Deitado -----

tERP_lying = squeeze(mean(targetEEG_lying, 1)); % ERP target (média do sinal na dimensão trials)
ntERP_lying = squeeze(mean(nontargetEEG_lying, 1)); % ERP não-target

%% Correção da linha de base

baseline = [-100;0];

baseidx = dsearchn(time',baseline);

% ----- Sentado -----

for e=1:length(channel_labels)
    tmn = mean(tERP_seated(baseidx(1):baseidx(2),e));   
    ntmn = mean(ntERP_seated(baseidx(1):baseidx(2),e));   

    tERP_seated(:,e) = tERP_seated(:,e) - tmn;
    ntERP_seated(:,e) = ntERP_seated(:,e) - ntmn;
end

% ----- fMRI -----

for e=1:length(channel_labels)
    tmn = mean(tERP_fmri(baseidx(1):baseidx(2),e));   
    ntmn = mean(ntERP_fmri(baseidx(1):baseidx(2),e));   

    tERP_fmri(:,e) = tERP_fmri(:,e) - tmn;
    ntERP_fmri(:,e) = ntERP_fmri(:,e) - ntmn;
end

% ----- Deitado -----

for e=1:length(channel_labels)
    tmn = mean(tERP_lying(baseidx(1):baseidx(2),e));   
    ntmn = mean(ntERP_lying(baseidx(1):baseidx(2),e));   

    tERP_lying(:,e) = tERP_lying(:,e) - tmn;
    ntERP_lying(:,e) = ntERP_lying(:,e) - ntmn;
end

%% Gráficos - Sentado

% Definição do tempo dos picos e respetivos índices 

tN170 = [95; 120]; % ms
idxN170 = dsearchn(time',tN170);

tN2 = [125; 160];
idxN2 = dsearchn(time',tN2);

tP3 = [250; 350];
idxP3 = dsearchn(time',tP3);

% % Gráficos dos ERPs no tempo
% 
% clim = [floor(min(min(tERP_seated))),ceil(max(max(tERP_seated)))];
% 
% aux = 0;
% 
% for i=1:length(channel_labels)
%     figure(fig)
%     hold on,
%     fig = fig + 1;
%     plot(time,tERP_seated(:,i),'b','LineWidth',2);
%     plot(time,ntERP_seated(:,i),'r','LineWidth',2);
%     xline(0, '--')    
%     ylim([(clim(1)+aux),(clim(2)-aux)])
%     xlim([min(time),max(time)])
%     title(channel_labels(i) + " - Sentado"), xlabel('Tempo (ms)'), ylabel('Amplitude (\muv)'), grid on
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
% 
%     hold off
% end
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
    plot(time,tERP_lying(:,i),'r','LineWidth',2);
    plot(time,tERP_fmri(:,i),'g','LineWidth',2);
    xline(0, '--')    
    ylim([(clim(1)+aux),(clim(2)-aux)])
    xlim([min(time),max(time)])
    title(channel_labels(i) + " - Target"), xlabel('Tempo (ms)'), ylabel('Amplitude (\muv)'), grid on
    
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

    legend('EEG Sentado','EEG Deitado','EEG c/ fMRI','Stimulus','Baseline',a, 'N2', 'P3','AutoUpdate', 'off')
    
    yline(0, ':')
    hold off
    
end

%% Comparação entre condições - Non-Target

% Gráficos dos ERPs no tempo

clim = [floor(min(min(ntERP_seated))),ceil(max(max(ntERP_seated))),
        floor(min(min(ntERP_fmri))),ceil(max(max(ntERP_fmri))),
        floor(min(min(ntERP_lying))),ceil(max(max(ntERP_lying)))];

clim = [min(min(clim)),max(max(clim))];

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

    legend('EEG Sentado','EEG Deitado','EEG c/ fMRI','Stimulus','Baseline',a, 'N2', 'P3','AutoUpdate', 'off')
    
    yline(0, ':')
    hold off
    
end

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

