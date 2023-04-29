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

data=readtable("S3_25_Male.csv");
data=table2array(data,1);

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

t = (-3:1/fs:3)';
s = 5*square(2*pi*t); % Sinal com a mesma frequência que o sinal da fotoresistência

sqwave = awgn(s, 20, 'measured'); % Gaussian white noise

plot(t, [sqwave, movmean(sqwave,10), movmedian(sqwave,20)]),
title('Verificar o efeito dos filtros na fase dos sinais'), 
xlabel('Tempo (s)'),
legend('Original','Média flutuante','Mediana flutuante');

% Resumindo, a mediana flutuante não distorce tanto as fases do sinal,
% principalmente no início das rising edges e falling edges. Isto é
% extremamente importante em análises como as dos ERPs pois é preciso
% saber com exatidão o instante em que o estímulo é apresentado.

%% Dados da fotoresistência 

photoVector=(data(:,end));

t = ((0:length(photoVector)-1)*1/fs)/60;
n_pnts = 20;
filt_photoVector = movmedian(photoVector,20);

% lpFilt = designfilt('lowpassfir','FilterOrder',1000,'PassbandFrequency',2, ...
%          'StopbandFrequency',2.5,'SampleRate',fs);
%fvtool(lpFilt)

%lpFilt_phVec = filtfilt(lpFilt,photoVector);

fig = 1;

% figure (fig); plot(t,photoVector), title('Fotoresistência'), xlabel('Tempo (min)');
% fig = fig + 1;
% 
% figure (fig); plot(t,filt_photoVector), title("Fotoresistência (Mediana flutuante com "+n_pnts+" pontos)"), xlabel('Tempo (min)');
% fig = fig + 1;

% figure (fig); plot(t,lpFilt_phVec), title("Fotoresistência (Filtro passa-baixo)"), xlabel('Tempo (min)');
% fig = fig + 1;

%% Re-referenciação do EEG 

eeg_data = data(:,1:end-1);

%car = bsxfun(@minus, eeg_data, mean(eeg_data,2));
%lap = laplacian_perrinX(eeg_data',[nch_locs.X],[nch_locs.Y],[nch_locs.Z]);

%% Pré-processamento do EEG

% Filtragem
% De acordo com o artigo com a literatura dos ERPs, um
% filtro passa-banda entre 0.5 e 30 Hz será usado.

% Band-pass
bpFilt = designfilt('bandpassfir','FilterOrder',1000, ...
         'CutoffFrequency1',0.5,'CutoffFrequency2',40, ...
         'SampleRate',fs); 
%fvtool(bpFilt) % Visualizar resposta em frequência e fase do filtro

%filt_data = filtfilt(bpFilt,car);
%filt_data = filtfilt(bpFilt,lap');
filt_data = filtfilt(bpFilt,eeg_data);

%% FFT do sinal 

% Verificar se o sinal é um EEG. Para isso, o sinal no domíno das
% frequências tem que seguir a forma 1/f.
% Isto após retirar as componentes DC e dos 50Hz da rede elétrica.

T = 1/fs; % Sampling period
L = length(filt_data(:,8)); % Length of signal
tm = (0:L-1)*T; % Time vector
Y = fft(filt_data(:,8)); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:L/2+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L; % Frequency vector

figure(fig)
fig = fig +1;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Retirar sub-harmónica dos 25 Hz

% Band-stop
bsFilt = designfilt('bandstopfir','FilterOrder',1000, ...
         'CutoffFrequency1',24.7,'CutoffFrequency2',25.3, ...
         'SampleRate',fs);
%fvtool(bsFilt)
filt_data = filtfilt(bsFilt,filt_data);

%% Deteção dos índices de segmentação

last_low = true; 
last_high = false;

falling_edges = [];
rising_edges = [];

for i = 1:length(filt_photoVector)
    photoresistor = filt_photoVector(i);
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

seg_target = [];
seg_nontarget = [];

for i = 1:length(falling_edges) % Verificar quais blocos pertecem aos targets e não targets
    if(min(filt_photoVector(falling_edges(i):rising_edges(i))) < 140)
        seg_nontarget = [seg_nontarget, falling_edges(i)];
    else
        seg_target = [seg_target, falling_edges(i)];
    end
end

seg_target = unique(seg_target);
seg_nontarget = unique(seg_nontarget);

%% Segmentação do EEG

bf = -0.100; % 100 ms antes do estímulo
af = 0.800; % 800 ms após o estímulo

idx_bf = floor(bf * fs);
idx_af = ceil(af * fs);

time = ((idx_bf:idx_af)/fs) * 1000; % Vetor do tempo

% Target

[targetEEG, nontargetEEG] = deal(zeros(length(seg_target), length(time), length(channel_labels))); % Pré-alocação dos dados

for i=1:length(seg_target)
    targetEEG(i,:,:) = cat(3,filt_data((seg_target(i)+idx_bf):(seg_target(i)+idx_af),:));
end

% Non-target

for i=1:length(seg_nontarget)
    nontargetEEG(i,:,:) = cat(3,filt_data((seg_nontarget(i)+idx_bf):(seg_nontarget(i)+idx_af),:));
end

%% Computação dos ERPs

tERP = squeeze(mean(targetEEG, 1)); % ERP target (média do sinal na dimensão trials)
ntERP = squeeze(mean(nontargetEEG, 1)); % ERP não-target

%% Correção da linha de base

baseline = [-100;0];

baseidx = dsearchn(time',baseline);

for e=1:length(channel_labels)
    tmn = mean(tERP(baseidx(1):baseidx(2),e));   
    ntmn = mean(ntERP(baseidx(1):baseidx(2),e));   

    tERP(:,e) = tERP(:,e) - tmn;
    ntERP(:,e) = ntERP(:,e) - ntmn;
end

%% Gráficos

% Definição do tempo dos picos e respetivos índices 

tN170 = [95; 120]; % ms
idxN170 = dsearchn(time',tN170);

tN2 = [125; 160];
idxN2 = dsearchn(time',tN2);

tP3 = [250; 350];
idxP3 = dsearchn(time',tP3);

% Gráficos dos ERPs no tempo

clim = [floor(min(min(tERP))),ceil(max(max(tERP)))];
aux = 0;

for i=1:length(channel_labels)
    figure(fig)
    hold on,
    fig = fig + 1;
    plot(time,tERP(:,i),'b','LineWidth',2);
    plot(time,ntERP(:,i),'r','LineWidth',2);
    xline(0, '--')    
    ylim([(clim(1)+aux),(clim(2)-aux)])
    xlim([min(time),max(time)])
    title(channel_labels(i)), xlabel('Tempo (ms)'), ylabel('Amplitude (\muv)'), grid on
    
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

    legend('Target','Non-target','Stimulus','Baseline',a, 'N2', 'P3')
   
    hold off
end

% Mapas topográficos

figure(fig)
fig = fig + 1;
tlo1=tiledlayout(1,3,'TileSpacing','tight','Padding','tight');

ax1=nexttile(tlo1);
topoplot(mean(tERP(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax1,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)

ax2=nexttile(tlo1);
topoplot(mean(tERP(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax2,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)

ax3=nexttile(tlo1);
topoplot(mean(tERP(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax3,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)

cbh1=colorbar(ax3);
cbh1.Layout.Tile = 'east';
cbh1.Label.String='\muV';
cbh1.Label.Rotation=0; 

figure(fig)
fig = fig + 1;
tlo2=tiledlayout(1,3,'TileSpacing','tight','Padding','tight');

ax4=nexttile(tlo2);
topoplot(mean(ntERP(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax4,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)

ax5=nexttile(tlo2);
topoplot(mean(ntERP(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax5,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)

ax6=nexttile(tlo2);
topoplot(mean(ntERP(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax6,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)

cbh1=colorbar(ax6);
cbh1.Layout.Tile = 'east';
cbh1.Label.String='\muV';
cbh1.Label.Rotation=0; 

