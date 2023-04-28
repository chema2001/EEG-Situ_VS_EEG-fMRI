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

%% Importar os dados
close all;
clc, clear;

addpath(genpath('eeglab2023.0')); % Path do EEGLAB
ch_locs=readlocs('BioSemi64.loc'); % Ficheiro com localização dos elétrodos

data=readtable("S1_21_Male.csv");
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

nch_locs = ch_locs(logical(idx));

%% Dados da fotoresistência 

photoVector=(data(:,end));
t = ((0:length(photoVector)-1)*1/fs)/60;

% Janela de tempo
window_size = 7;

% Padding do sinal
padded_signal = [zeros(window_size-1, 1); photoVector; zeros(window_size-1, 1)];
filt_photoVector = zeros(size(photoVector));

for i = 1:length(photoVector)
    window = padded_signal(i:i+window_size-1);
    
    % Valor da mediana
    median_value = median(window);
    filt_photoVector(i) = median_value;   
end

fig = 1;
figure (fig); plot(t,photoVector), title('Fotoresistência'), xlabel('Tempo (min)');
fig = fig + 1;

figure (fig); plot(t,filt_photoVector), title('Fotoresistência (filtrado)'), xlabel('Tempo (min)');
fig = fig + 1;

%% Pré-processamento do EEG
% Filtragem
% De acordo com o artigo 'A Rapid Face Recognition BCI System using
% Single-trial ERP', um filtro passa-banda entre 0.5 e 30 Hz será usado

% Band-pass
bpFilt = designfilt('bandpassfir','FilterOrder',1000, ...
         'CutoffFrequency1',0.5,'CutoffFrequency2',30, ...
         'SampleRate',fs); 
%fvtool(bpFilt) % Visualizar resposta em frequência e fase do filtro
filt_data = filtfilt(bpFilt,data(:,1:end-1));

% Band-stop
bsFilt = designfilt('bandstopfir','FilterOrder',1000, ...
         'CutoffFrequency1',24.5,'CutoffFrequency2',25.5, ...
         'SampleRate',fs);
%fvtool(bsFilt)
filt_data = filtfilt(bsFilt,filt_data);

%% Indices da segmentação do EEG target
% Sinal da fotoresistência não filtrado 

th=230;
targetVector=double((photoVector>th));

% locs -> Índice da primeira amostra
% dur -> duração dos blocos em nº de amostras
[~,locs,dur,~] = findpeaks(targetVector);

figure(fig), histogram((dur/fs)*1000,60), title('Duração dos blocos target'), xlabel('Tempo (s)') % Verificar o tamanho dos blocos (trials)
fig = fig + 1;

idx = ((dur/fs)*1000 <= 480); % Retira os blocos com nº de amostras/duração inválida (<470 ms)
dur(idx) = [];
locs(idx) = [];
minw = min(dur);
dur(1:end) = minw;

% Sinal da fotoresistência filtrado 

targetVector=double((filt_photoVector>th));

% locs -> Índice da primeira amostra
% dur -> duração dos blocos em nº de amostras
[~,f_locs,f_dur,~] = findpeaks(targetVector);

figure(fig), histogram((f_dur/fs)*1000,60), title('Duração dos blocos target'), xlabel('Tempo (s)') % Verificar o tamanho dos blocos (trials)
fig = fig + 1;

f_idx = ((f_dur/fs)*1000 <= 480); % Retira os blocos com nº de amostras/duração inválida (<470 ms)
f_dur(f_idx) = [];
f_locs(f_idx) = [];
f_minw = min(f_dur);
f_dur(1:end) = f_minw;

shift = abs(locs - f_locs);

f_locs = f_locs - shift; % Ajustar a mudança de fase que ocorre com o filtro de mediana, o qual é constante para



%% Segmentação do EEG target

targetEEG = [];

for i=1:length(dur)
    targetEEG = [targetEEG; filt_data(locs(i)-1:locs(i)+f_dur(i),:)];
end

segmented_targetEEG = reshape(targetEEG, f_dur(1)+2, length(f_dur),length(channel_labels));

%% Indices da segmentação do EEG não-target

th = 160;
non_targetVector=double((filt_photoVector>th) & ~targetVector);

% locs -> Índice da primeira amostra
% dur -> duração dos blocos em nº de amostras
[~,nont_locs,nont_dur,~] = findpeaks(non_targetVector);

figure(fig), histogram((nont_dur/fs)*1000,60), title('Duração dos blocos non-target'), xlabel('Tempo (s)') % Verificar o tamanho dos blocos (trials)
fig = fig + 1;

idx = ((nont_dur/fs)*1000 <= 480); % Retira os blocos com nº de amostras/duração inválida (<470 ms)
nont_dur(idx) = [];
nont_locs(idx) = [];
minw = min(nont_dur);
nont_dur(1:end) = minw;

nont_locs = nont_locs - shift(1);

%% Segmentação do EEG target

nontargetEEG = [];

for i=1:length(nont_dur)
    nontargetEEG = [nontargetEEG; filt_data(nont_locs(i)-1:nont_locs(i)+nont_dur(i),:)];
end

segmented_nontargetEEG = reshape(nontargetEEG, nont_dur(1)+2, length(nont_dur),length(channel_labels));

%% Computação dos ERPs

tg_erp = squeeze(mean(segmented_targetEEG, 2)); % Obtenção do ERP target (média do sinal na dimensão trials)
ntg_erp = squeeze(mean(segmented_nontargetEEG, 2)); % Obtenção do ERP não-target

t = (linspace(0,(dur(1)+2)/fs,dur(1)+2)).* 1000;

% Definição do tempo dos picos e respetivos índices 

tN170 = [95; 135]; % ms
idxN170 = dsearchn(t',tN170);

tN2 = [180; 230];
idxN2 = dsearchn(t',tN2);

tP3 = [325; 450];
idxP3 = dsearchn(t',tP3);

% Gráficos dos ERPs no tempo

clim = [floor(min(min(tg_erp))),ceil(max(max(tg_erp)))];
aux = 0;

for i=1:length(channel_labels)
    figure(fig)
    hold on,
    fig = fig + 1;
    plot(t,tg_erp(:,i),'b','LineWidth',2);
    plot(t,ntg_erp(:,i),'r','LineWidth',2);
    ylim([(clim(1)+aux),(clim(2)-aux)])
    xlim([min(t),max(t)])
    title(channel_labels(i)), xlabel('Tempo (ms)'), ylabel('Amplitude (\muv)'), yline(0, '--'), grid on
    
    if(contains(channel_labels(i),"z") || (contains(channel_labels(i),"C") && ~(contains(channel_labels(i),"P"))) ) % VPP
        c = 'blue';
        a = 'VPP';
        N2 = [tN2(1)-30, tN2(2)-30, tN2(2)-10, tN2(1)-10];
    else
        c = 'cyan';
        a = 'N170';
        N2 = [tN2(1), tN2(2), tN2(2), tN2(1)];
    end
   
    N170 = [tN170(1), tN170(2), tN170(2), tN170(1)];
    P3 = [tP3(1), tP3(2), tP3(2), tP3(1)];

    % Estas patches precisam de ser ajustadas

    ph1 = patch(N170,[clim(1) clim(1) clim(2) clim(2)],c);
    set(ph1,'facealpha',.1,'edgecolor','none')

    ph2 = patch(N2,[clim(1) clim(1) clim(2) clim(2)],'y');
    set(ph2,'facealpha',.1,'edgecolor','none')
    
    ph3 = patch(P3,[clim(1) clim(1) clim(2) clim(2)],'r');
    set(ph3,'facealpha',.1,'edgecolor','none')

    legend('Target','Non-target','Baseline','N170 - VPP', 'N2', 'P3')
   
    hold off
end

% Mapas topográficos

figure(fig)
fig = fig + 1;
tlo1=tiledlayout(1,3,'TileSpacing','tight','Padding','tight');

ax1=nexttile(tlo1);
topoplot(mean(tg_erp(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax1,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)

ax2=nexttile(tlo1);
topoplot(mean(tg_erp(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax2,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)

ax3=nexttile(tlo1);
topoplot(mean(tg_erp(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax3,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)

cbh1=colorbar(ax3);
cbh1.Layout.Tile = 'east';
cbh1.Label.String='\muV';
cbh1.Label.Rotation=0; 

figure(fig)
fig = fig + 1;
tlo2=tiledlayout(1,3,'TileSpacing','tight','Padding','tight');

ax4=nexttile(tlo2);
topoplot(mean(ntg_erp(idxN170(1):idxN170(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax4,"N170 - VPP (" + tN170(1) +" - "+ tN170(2) + " ms)"  ,'FontSize',12)

ax5=nexttile(tlo2);
topoplot(mean(ntg_erp(idxN2(1):idxN2(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax5,"N2 (" + tN2(1) +" - "+ tN2(2) + " ms)",'FontSize',12)

ax6=nexttile(tlo2);
topoplot(mean(ntg_erp(idxP3(1):idxP3(2),:)),nch_locs,'maplimits',clim,'numcontour',6,'electrodes','on');
title(ax6,"P3 (" + tP3(1) +" - "+ tP3(2) + " ms)",'FontSize',12)

cbh1=colorbar(ax6);
cbh1.Layout.Tile = 'east';
cbh1.Label.String='\muV';
cbh1.Label.Rotation=0; 