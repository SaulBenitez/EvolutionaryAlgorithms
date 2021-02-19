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
%--------ED
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

%--------CPSO
C1 = 2.05;
C2 = 2.05;
w = 0.729;
vmin = -0.1*xmax;
vmax = 0.1*xmax;      
p = 0.7;
rho = 0.1;

%% Funcion objetivo
f = @Ackley;

%% Ejecucion de los algoritmos
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
    
    fprintf('************************* Corrida %d - CPSO *************************\n', run);
    CPSO(f, NP, D, DP, G, NEJ, xmin, xmax, C1, C2, w, vmin, vmax, p, rho, fprint, run);
end

%% Guardado de parametros
save(strcat('All','_',func2str(f),'_Config'))
