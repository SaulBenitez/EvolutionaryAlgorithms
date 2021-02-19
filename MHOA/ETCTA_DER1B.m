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

%% Parametros de impresion resultados
fprint = 1;
 
%% Bound design variables
%-------------------------------- varrho limits --------------------------------
varrhoMax = [90.0e3, 30.0e3, 30e3, 35, 30, 10];     varrhoMin = [0, 0, 0, 0, 0, 0];
%-------------------------------- Epsilon limits --------------------------------
epsMax = [4, 4, 5];                                 epsMin = [0, 0, 0];
%-------------------------------- Delta limits --------------------------------
DeltaMax = [3, 2, 3];                               DeltaMin = [0, 0, 0];
%-------------------------------- sigma limits --------------------------------
sigmaMax = 1;                                       sigmaMin = 0;
%-------------------------------- omega limits --------------------------------
omegaMax = 2;                                       omegaMin = 0;
%-------------------------------- psi limits --------------------------------
psiMax = [8, 6, 2];                                 psiMin = [0, 0, 0];

%% Parametros del problema
NEI = 30;
NA = 12;
D = 17;
DP = 4;
NP = 50;
G = 2000;
NEJ = NP*G;
xmin = [varrhoMin, epsMin, DeltaMin, sigmaMin, omegaMin, psiMin];
xmax = [varrhoMax, epsMax, DeltaMax, sigmaMax, omegaMax, psiMax];
xi = [1000, 0.01];

%% Parametros del algoritmo
Fmin = 0.1;
Fmax = 0.9;
Kmin = 0.1;
Kmax = 0.9;
CR = 0.6;

%% Ejecucion de los algoritmos
f = @ETCTA;
parfor run = 1:NEI
    fprintf('************************* Corrida %d - Rand/1/Bin *************************\n', run);
    DER1B(f, NP, G, D, DP, Fmin, Fmax, CR, NEJ, xmin, xmax, fprint, run)
end

%% Guardado de parametros
save(strcat('DER1B','_',func2str(f),'_Config'))
