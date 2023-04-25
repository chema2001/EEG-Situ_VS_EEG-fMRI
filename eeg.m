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

addpath(genpath('C:\Users\migue\OneDrive\Ambiente de Trabalho\EEG stuff\eeglab2023.0')); % Path do EEGLAB
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

fig = 1;
figure (fig); plot(t,photoVector), title('Fotoresistência'), xlabel('Tempo (min)');
fig = fig + 1;

%% Pré-processamento
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

%% Segmentação do EEG target

th=230;
targetVector=double((photoVector>th));

% locs -> Índice da primeira amostra
% dur -> duração dos blocos em nº de amostras
[~,locs,dur,~] = findpeaks(targetVector); 

figure(fig), histogram((dur/fs)*1000,60), title('Duração dos blocos target'), xlabel('Tempo (s)') % Verificar o tamanho dos blocos (trials)
fig = fig + 1;

idx = ((dur/fs)*1000 <= 470); % Retira os blocos com nº de amostras/duração inválida (<470 ms)
dur(idx) = [];
minw = min(dur);
dur(1:end) = minw;

targetEEG = [];

for i=1:length(dur)
    targetEEG = [targetEEG; filt_data(locs(i)-1:locs(i)+dur(i),:)];
end

segmented_targetEEG = reshape(targetEEG, dur(1)+2, length(dur),length(channel_labels));

tg_erp = squeeze(mean(segmented_targetEEG, 2)); % Obtenção do ERP (média do sinal na dimensão trials)

t = (linspace(0,(dur(1)+2)/fs,dur(1)+2)).* 1000;

% Definição do tempo dos picos e respetivos índices 

tN170 = [90; 145]; % ms
idxN170 = dsearchn(t',tN170);

tN2 = [180; 230];
idxN2 = dsearchn(t',tN2);

tP3 = [325; 450];
idxP3 = dsearchn(t',tP3);

% Gráficos dos ERPs no tempo

clim = [floor(min(min(tg_erp))),ceil(max(max(tg_erp)))];
aux = 3;

for i=1:length(channel_labels)
    figure(fig)
    hold on,
    fig = fig + 1;
    plot(t,tg_erp(:,i),'LineWidth',2); 
    ylim([(clim(1)+aux),(clim(2)-aux)])
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

%     [~,h_legend] = legend(a,'N2','P3');
%     PatchInLegend = findobj(h_legend, 'type', 'patch');
%     set(PatchInLegend(1));
%     set(PatchInLegend(2));
%     set(PatchInLegend(3));
   
    hold off
end

% Mapas topográficos

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

%% Segmentação do EEG target não target

th1 = 170;
nontargetVector=double(photoVector>=th1 & ~targetVector);

[~,locs,w,~] = findpeaks(nontargetVector);

%%
% 
% targetEEGs={};
% channel16Vector=data(:,16);
% for i=1:length(locs)
%     targetEEGs{i}=channel16Vector(locs(i):locs(i)+w(i)).';
% end
% 
% for i=1:length(targetEEGs)
%     targetEEGs{i}=[targetEEGs{i},nan(1,max(w)-length(targetEEGs{i}))];
% end
% 
% b=cell2mat(targetEEGs);
% c = mean(b,1,'omitnan');
% ind = ~isnan(c);
% c = c(ind);
% c = highpass(c, 0.5, 125);
% c = lowpass(c, 35, 125);
% 
% d=linspace(0,100,length(c));
% figure(3); plot(d,c);
