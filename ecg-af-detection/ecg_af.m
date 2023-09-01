%% 1. DETECCIÓN DE FIBRILACIÓN AURICULAR MEDIANTE CLUSTERING 
clear; close all; clc;

%Zero-Padding for Border Effects
dwtmode('zpd','nodisp');

%Cargamos nuestra señal (1 = normal, 2 = AF)
load ecg;

%% 1.1 Pre-procesamiento

x_ax = [1000 1005];           % Xlim(s)

Fs = 256;                     % Frecuencia de muestreo
t = (1:1:length(ecg))*(1/Fs); % Tiempo de muestreo

%% 1.2. Análisis en el dominio de la frecuencia y filtros

F = abs(fft(ecg));  % Transformada de Fourier y abs value
F = F(1:end/2+1);   % Analizamos la 1ra mitad (fft es simétrica)
F = (F/max(F));     % Normalizamos los valores
L = length(ecg);    % Longitud de la señal
f = Fs*(0:L/2)/L;   % xlim(f)


% 1.2.1. Valores cercanos a cero

% HIGH PASS
N     = 20;   % Order
Fstop = 0.1;  % Stopband Frequency
Astop = 80;   % Stopband Attenuation (dB)

% Construcción del filtro IRR
h1  = fdesign.highpass('N,Fst,Ast', N, Fstop, Astop, Fs);
Hd1 = design(h1, 'cheby2');

% Señal luego del 1er filtro
ecg_filt1 = filter(Hd1, ecg);

% NOTCH/BANDSTOP

N      = 20;  % Order
Fstop1 = 49;  % First Stopband Frequency
Fstop2 = 51;  % Second Stopband Frequency
Astop  = 100; % Stopband Attenuation (dB)

% Construcción del filtro IRR
h2  = fdesign.bandstop('N,Fst1,Fst2,Ast', N, Fstop1, Fstop2, Astop, Fs);
Hd2 = design(h2, 'cheby2');

% Señal luego del 2do filtro
ecg_filt2 = filter(Hd2, ecg_filt1);


% MOVING AVERAGE FILTER
% Filtro para suavizar la señal
ord = 30;                            % Order
b = (1/ord)*(ones(1,ord));           % Filter 
ECG_trend = filter(b,[1],ecg_filt2); % Filter apply
ECG_clean = ecg_filt2-ECG_trend;     % Clean signal

% Análisis en el dominio de la frecuencia (Señal limpia).

F1 = abs(fft(ECG_clean));  % Transformada de Fourier y abs value
F1 = F1(1:end/2+1);        % Analizamos la 1ra mitad (fft es simétrica)
F1 = (F1/max(F1));         % Normalizamos los valores
L1 = length(ecg);          % Longitud de la señal
f1 = Fs*(0:L1/2)/L1;       % xlim(f)

% Graficamos

figure;

subplot(2,2,1);
plot(t,ecg);
xlim(x_ax);
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');
title('ECG original');

subplot(2,2,2);
plot(f,F);
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
title('Transformada rápida de Fourier - ECG original');

subplot(2,2,3);
plot(t,ECG_clean);
xlim(x_ax);
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');
title('ECG filtrada');

subplot(2,2,4);
plot(f1,F1);
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
title('Transformada rápida de Fourier - ECG filtrada')

%% 1.3. Segmentación

s10 = Fs*10;              % Muestras cada 10s
n = int64(L/s10)-1;       % Aproximación del número de segmentos

% Matriz (rows = grupo de 10s, cols = muestras)
ecg_seg = ones(n,s10);
t_10 = (1:1:s10)*(1/Fs);  % Vector de tiempo para 10s
label = []; 
 
for i = 1:n
    ecg_seg(i,:) = ECG_clean(i*s10+1:(i+1)*s10);                 % Llenamos la matriz con la señal
    label = [label int64(mean(annotations(i*s10+1:(i+1)*s10)))]; % Concatenamos valores
end

% Graficamos
figure;

% Visualización de una señal ECG normal
subplot(2,1,1);
plot(t_10,ecg_seg(150,:))
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');
title('ECG Normal');

% Visualización de AF, debido a la variación en RR distance
subplot(2,1,2);
plot(t_10,ecg_seg(20,:))
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');
title('ECG AF');

%% 1.4 Extracción y visualización de características

features = ones(n,17); % Matriz de características
not_good = [];         % Guardamos el índice de las anomalías


for j = 1:n
    
    tmp_ecg = ecg_seg(j,:);
    
    % Energy -----------------------------------
    
    Eng = sum(tmp_ecg.^2);
    
    % Wavelet Transform----------------------------------

    wname = 'db4';                    % Daubechies wavelet order four
    [c,l] = wavedec(tmp_ecg,5,wname); % Wavelet decomposition 5 levels
    approx = appcoef(c,l,'db4');      % Approximation coefficients
    [cd1,cd2,cd3,cd4,cd5] = detcoef(c,l,[1 2 3 4 5]);
    
    % Reconstruct wavelets
    a5_s = wrcoef('a',c,l,'db4');
    cd1_s = wrcoef('d',c,l,'db4',1);
    cd2_s = wrcoef('d',c,l,'db4',2);
    cd3_s = wrcoef('d',c,l,'db4',3);
    cd4_s = wrcoef('d',c,l,'db4',4);
    cd5_s = wrcoef('d',c,l,'db4',5);

    % Absolute wavelet subband energies
    E1 = [];
    E2 = [];
    E3 = [];
    E4 = [];
    E5 = [];
    E0 = [];
    
    for i=1 : length(a5_s)
        E0(end+1)=(a5_s(i)^2);
    end
    E0=sum(E0);
    
    for i=1 : length(cd1_s)
        E1(end+1)=(cd1_s(i)^2);
    end
    E1=sum(E1);
    
    for i=1 : length(cd2_s)
        E2(end+1)=(cd2_s(i)^2);
    end
    E2=sum(E2);
    
    for i=1 : length(cd3_s)
        E3(end+1)=(cd3_s(i)^2);
    end
    E3=sum(E3);
    
    for i=1 : length(cd4_s)
        E4(end+1)=(cd4_s(i)^2);
    end
    E4=sum(E4);
    
    for i=1 : length(cd5_s)
        E5(end+1)=(cd5_s(i)^2);
    end
    E5=sum(E5);
    
    % Total energy from levels
    Ewavelet = E0 + E1 + E2 + E3 + E4 + E5;

    % Relative wavelet subband energies
    ER0 = E0/Ewavelet;
    ER1 = E1/Ewavelet;
    ER2 = E2/Ewavelet;
    ER3 = E3/Ewavelet;
    ER4 = E4/Ewavelet;
    ER5 = E5/Ewavelet;

    % HRV  -------------------------------------  

    % QRS detection
    
    tmp_ecg_norm = tmp_ecg*(1/max(tmp_ecg)); % Normalizamos la data
    tmp_norm_2 = (tmp_ecg_norm+2).^2;        % Definimos los picos a detectar
    [pks, locs]= findpeaks(tmp_norm_2,Fs,'MinPeakHeight',1.8*mean(tmp_norm_2),'MinPeakWidth',(1/Fs)*5); % Encontramos picos en la señal
    
    % HRV and Interpolation
    
    HRV = diff(locs); % Encontramos la HRV de la señal
    pks1 = [];
    for i = 1:length(locs)
        pks1 = [pks1 tmp_ecg(locs(i)*Fs+1)];
    end
    
    t1 = (1:1:length(HRV))*(10/length(HRV));
    
    % RMSSD -------------------------------------
    
    sum1 = 0;
    for i = 1:length(HRV)-1
        sum1 = sum1 + (HRV(i+1)-HRV(i))^2;
    end   
    RMSSD = sqrt((1/(length(HRV)))*sum1);
    
    % pRR50 -------------------------------------
    
    sum2 = 0;
    for i = 1:length(HRV)-1
        dRRi = abs(HRV(i+1)-HRV(i));
        if dRRi > 0.5
            sum2 = sum2 + 1;
        end
    end   
    
    pRR50 = sum2/(length(HRV)-1);
  
    % minRR -------------------------------------
    
    minRR = min(HRV);
    
    
    % GRAFICAMOS ---------------------
   
    if j == 150 || j == 20
        
        figure;
        
        subplot(2,1,1);
        hold;
        plot(t_10,tmp_ecg);
        scatter(locs,pks1);
        hold off;
        
        xlabel('Tiempo (s)');
        ylabel('Amplitud (mV)');
        
        if j == 150
            title('QRS Detection (NORMAL)');
        else
            title('QRS Detection (AF)');
        end

        subplot(2,1,2);
        plot(t1, HRV);
        xlabel('Tiempo (s)');
        ylabel('RRI (s)');
        
        if j == 150
            title('Tacograma (NORMAL)');
        else
            title('Tacograma (AF)');
        end
        
    end
   
    % Llenamos nuestra matriz de características
    if isnan(RMSSD)
        not_good = [not_good j];
    else
        features(j,:) = [E0,E1,E2,E3,E4,E5,...
                         ER0,ER1,ER2,ER3,ER4,ER5...
                         Ewavelet,Eng,RMSSD,pRR50,minRR];
    end
    
end

clc;

% Corrección de segmentos
features(not_good,:) = [];

% Normalización
norm_feat = features./max(features);

feat_mean = mean(features); % Get mean per each column
feat_sd = std(features);    % Get std per column

v_names = {'E0','E1','E2','E3','E4','E5','ER0','ER1','ER2','ER3','ER4','ER5'...
                         ,'Ewavelet','Energy','RMSSD','pRR50','minRR'};
T = table(feat_mean',feat_sd','VariableNames', {'Mean','Standard_Deviation'},...
    'RowName', v_names);

disp(T);

% Reducción de las dimensiones
cov_mat = cov(norm_feat);
eig_val = fliplr(eig(cov_mat)'); %fliplr for arragne vect from max to min

figure;
plot(eig_val);
xlabel('Components');
ylabel('Eigenvalues');
title('Eigenvalues of Features')

% PCA
[coef, score] = pca(norm_feat,'NumComponents',2);

% Trazado de características
cmp = 2;

figure;
biplot(coef(:,1:cmp),'Scores',score(:,1:cmp),'VarLabels',v_names)

[coef1, score1, ] = pca(norm_feat(:,10:17),'NumComponents',2);

figure;
biplot(coef1(:,1:cmp),'Scores',score1(:,1:cmp),'VarLabels',{v_names{10:17}})


%% 1.5 Clasificación no supervisada - Clustering

% K means
sc1 = score1(:,1:2);

% Elegimos el mejor cluster de 5 oportunidades
[idx,C] = kmeans(sc1,2,'Replicates',5);

% 1.5.1 Evaluación del modelo
 
% Precisión
n = length(idx);

TP = 0;
FP = 0;
TN = 0;
FN = 0;

% 1 = normal ECG
% 2 = AF 

% Calculamos las métricas
for i = 1:n
    if (idx(i) == 2 & label(i) == 2)
      TP = TP + 1;
    elseif (idx(i) == 1 & label(i) == 2)
      FN = FN + 1;
    elseif (idx(i) == 1 & label(i) == 1)
      TN = TN + 1;
    else
      FP = FP + 1;
    end
end


acc = (TP + TN)/n;
sens = (TP/(TP + FN));
spec = (TN/(TN + FP));

% Evitamos posibles errores del algoritmo no supervisado
if (acc < 0.5)
    idx = 3 - idx;
end

TP = 0;
FP = 0;
TN = 0;
FN = 0;

for i = 1:n
    if (idx(i) == 2 & label(i) == 2)
      TP = TP + 1;
    elseif (idx(i) == 1 & label(i) == 2)
      FN = FN + 1;
    elseif (idx(i) == 1 & label(i) == 1)
      TN = TN + 1;
    else
      FP = FP + 1;
    end
end


acc = (TP + TN)/n;
sens = (TP/(TP + FN));
spec = (TN/(TN + FP));



T1 = table(acc,sens,spec,'VariableNames', {'Accuracy','Sensitivfivity','Specificity'},...
    'RowName', {'values'});
disp(T1);

% Graficamos
figure;

h = biplot(coef1(:,1:cmp),'Scores',score1(:,1:cmp),'VarLabels',{v_names{10:17}});
h_type = get(h, 'tag');
dot = h(strcmp(h_type,'obsmarker'));  

for i = 1:length(dot)
    if idx(i) == 1
        set(dot(i),'Color','green')
    else
        set(dot(i),'Color','red')
    end
end

disp('done')