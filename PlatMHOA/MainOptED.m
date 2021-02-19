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
G = 5;
NEJ = NP*G;
xmin = [varrhoMin, epsMin, DeltaMin, sigmaMin, omegaMin, psiMin];
xmax = [varrhoMax, epsMax, DeltaMax, sigmaMax, omegaMax, psiMax];
xi = [1000, 0.01];

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

%% Ejecucion de los algoritmos
f = @ETCTA;
parfor run = 1:NEI
    fprintf('************************* Corrida %d - Best/1/Bin *************************\n', run);
    DEB1B(f, NP, G, D, DP, Fmin, Fmax, CR, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Best/1/Exp *************************\n', run);
    DEB1E(f, NP, G, D, DP, Fmin, Fmax, CR, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Current-to-Best/1 *************************\n', run);
    DECB1(f, NP, G, D, DP, Fmin, Fmax, Kmin, Kmax, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Current-to-Best/1/Bin *************************\n', run);
    DECB1B(f, NP, G, D, DP, Fmin, Fmax, Kmin, Kmax, CR, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Current-to-Best/1/Exp *************************\n', run);
    DECB1E(f, NP, G, D, DP, Fmin, Fmax, Kmin, Kmax, CR, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Current-to-Rand/1 *************************\n', run);
    DECR1(f, NP, G, D, DP, Fmin, Fmax, Kmin, Kmax, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Current-to-Rand/1/Bin *************************\n', run);
    DECR1B(f, NP, G, D, DP, Fmin, Fmax, Kmin, Kmax, CR, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Rand/1/Exp *************************\n', run);
    DECR1E(f, NP, G, D, DP, Fmin, Fmax, Kmin, Kmax, CR, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Rand/1/Bin *************************\n', run);
    DER1B(f, NP, G, D, DP, Fmin, Fmax, CR, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Current-to-Rand/1/Exp *************************\n', run);
    DER1E(f, NP, G, D, DP, Fmin, Fmax, CR, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Rand/1/Either-or *************************\n', run);
    DER1EO(f, NP, G, D, DP, Fmin, Fmax, Kmin, Kmax, PF, NEJ, xmin, xmax, fprint, run)

    fprintf('************************* Corrida %d - Rand/2/Dir *************************\n', run);
    DER2D(f, NP, G, D, DP, Fmin, Fmax, NEJ, xmin, xmax, fprint, run)
end

%% Guardado de parametros
save(strcat('All','_',func2str(f),'_Config'))
