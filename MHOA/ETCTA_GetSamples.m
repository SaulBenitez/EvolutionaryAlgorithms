%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autor:    Saul Enrique Benitez Garcia                     %
% Archivo:  ParametricStatistics                            %
% Titulo:   Gather data sample from algorithms results      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%********************************************************************
% Description: Gather the best individuals for statistics
%   - POP: Best population of all independent execution
%   - POPhb: History of the best individual for each algorithm
%********************************************************************

clear; close all; format short; clc;

addpath('Statistics')

%% Parametros del problema
NEI = 30;
D = 17;
DP = 4;
NP = 50;
G = 2000;
NEJ = NP*G;
NAE = 14; %Numero de algoritmos evaluados
DMANAME = ["DER1B";
            "DER1E"; 
            "DEB1B"; 
            "DEB1E"; 
            "DECR1"; 
            "DECR1B"; 
            "DECR1E"; 
            "DECB1"; 
            "DECB1B"; 
            "DECB1E"; 
            "DER2D"; 
            "DER1EO"; 
            "GA"; 
            "CPSO"];
MANAME = ["DE/Rand/1/Bin"; 
            "DE/Rand/1/Exp"; 
            "DE/Best/1/Bin"; 
            "DE/Best/1/Exp"; 
            "DE/Current-to-Rand/1"; 
            "DE/Current-to-Rand/1/Bin"; 
            "DE/Current-to-Rand/1/Exp"; 
            "DE/Current-to-Best/1"; 
            "DE/Current-to-Best/1/Bin"; 
            "DE/Current-to-Best/1/Exp"; 
            "DE/Rand/2/Dir"; 
            "DE/Rand/1/Either-or";
            "GA";
            "CPSO"];
f = @ETCTA;
FUNCNAME = func2str(f);

%% Load data of files for each algorith
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datrfldr = 'Data'; % Root data folder
numfldr =  1; % Number folder where data execution will be evaluated +-
for i = 1:NAE
    algfldrs{i} = strcat(DMANAME(i),'_',num2str(numfldr));
    datapath{i} = fullfile(datrfldr,algfldrs{i});
    dataname(i) = strcat(DMANAME(i),'_',FUNCNAME,'_D',num2str(D),'_');
    datafullpath{i} = fullfile(datapath{i},dataname(i));
end

file_names = {'DE_Rand_1_Bin'; 
            'DE_Rand_1_Exp'; 
            'DE_Best_1_Bin'; 
            'DE_Best_1_Exp'; 
            'DE_Current-to-Rand_1'; 
            'DE_Current-to-Rand_1_Bin'; 
            'DE_Current-to-Rand_1_Exp'; 
            'DE_Current-to-Best_1'; 
            'DE_Current-to-Best_1_Bin'; 
            'DE_Current-to-Best_1_Exp'; 
            'DE_Rand_2_Dir'; 
            'DE_Rand_1_Either-or';
            'GA'
            'CPSO'
            };


%% Almacenamiento de datos de la última generacion de cada algoritmo evaluado
%------------ Extraccion de datos de variantes de ED -------------       
XP = zeros(NP,D+DP,NEI,NAE);
XPhb = zeros(G,D+DP,NAE);
for i = 1:NAE
    for j = 1:NEI
        load(strcat(datafullpath{i},num2str(j),'.mat'), 'AP')
        XP(:,:,j,i) = AP(:,:,end);
    end
end
for i = 1:NAE
    [~, IDXmin] = min(XP(:,D+1,:,i));
    [~, IDXminexe] = min(min(XP(:,D+1,:,i)));
    load(strcat(datafullpath{i},num2str(IDXminexe),'.mat'), 'AP')
    for k = 1:G
        XPhb(k,:,i) = AP(IDXmin(IDXminexe),:,k);
    end
end

POP = XP;
POPhb = XPhb;
clear AP XP XPhb DMANAME FUNCNAME MANAME
clear IDXmin IDXminexe i j k
clear algfldrs datafullpath dataname datapath datrfldr numfldr


save('Statistics/ETC-AlgorithmsDataSamples')
