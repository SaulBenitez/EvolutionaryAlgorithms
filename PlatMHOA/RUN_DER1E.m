%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autor:    Saul Enrique Benitez Garcia                     %
% Archivo:  MainOpt.m                                       %
% Titulo:   Evolucion diferencial 12 variantes              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%********************************************************************
% Descripcion: 
%********************************************************************
clear; close all; format long; clc;

% rng('shuffle');
% warning('off');

%% Additional paths
addpath('Algorithms')
addpath('Problems')

%% Parametros de impresion resultados
fprint = 1;

%% Parametros del problema
NEI = 1;
NA = 12;
D = 20;
DP = 4;
NP = 50;
G = 1000;
NEJ = NP*G;
xmin = -15*ones(1,D);
xmax = 30*ones(1,D);

%% Parametros del algoritmo
F = [0.3,0.3,0.9];
Fmin = 0.1;
Fmax = 0.9;
Kmin = 0.1;
Kmax = 0.9;
K = [0.3, 0.3, 0.9];
PF = 0.5;
CR = 0.6;
CRmin = 0.1;
CRmax = 0.9;

%% Funcion objetivo
f = @Ackley;

%% Ejecucion del algoritmo
parfor run = 1:NEI
    fprintf('************************* Corrida %d - Rand/1/Exp *************************\n', run);
    DER1E(f, NP, G, D, DP, Fmin, Fmax, CR, NEJ, xmin, xmax, fprint, run)
end

%% Guardado de parametros
save(strcat('DER1E','_',func2str(f),'_Config'))
