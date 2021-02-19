%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autor:    Saul Enrique Benitez Garcia                     %
% Archivo:  pso.m                                           %
% Titulo:   Particle Swarm Optimization                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%****************************************************************
% Variables                                                     *
%   NC -> Numero de corridas                                    *
%   D -> Numero de variables de diseño                          *
%   G -> Numero de generaciones                                 *
%   NP -> Numero de individuos                                  *
%   F -> Factor de escala                                       *
%   CR -> Factor de cruza                                       *
%   Fmin -> Factor de escala minimo                             *
%   Fmax -> Factor de escala maximo                             *
%   NEJ -> Numero de evaluaciones maximo de J                   *
%   J -> valor de funcion objetivo                              *
%   Jprom -> Promedio de J por cada generacion                  *
%   Jstd -> Desciacion estandar de J por cada generacion        *
%   Jbest -> Mejor J                                            *
%   Jworst -> Peor J                                            *
%   x -> Poblacion                                              *
%   Gk -> Generacion en que CEJ llego a NEJ                     *
%   xmin -> valor minimo de x                                   *
%   xmin -> valor maximo de x                                   *
%****************************************************************

clear; close all; clc

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
D = 17;
DP = 4;
NP = 50;
G = 2000;
NEJ = NP*G;
xmin = [varrhoMin, epsMin, DeltaMin, sigmaMin, omegaMin, psiMin];
xmax = [varrhoMax, epsMax, DeltaMax, sigmaMax, omegaMax, psiMax];
xi = [1000, 0.01];

%% Parametros del algoritmo
C1 = 2.05;
C2 = 2.05;
w = 0.729;
vmin = -0.1*xmax;
vmax = 0.1*xmax;      
p = 0.7;
rho = 0.1;

%% Funcion objetivo
f = @ETCTA;

%% Ejecucion de algoritmo
parfor run = 1:NEI
    fprintf('************************* Corrida %d - CPSO *************************\n', run);
    CPSO(f, NP, D, DP, G, NEJ, xmin, xmax, C1, C2, w, vmin, vmax, p, rho, fprint, run);
end

%% Guardado de parametros
save(strcat('CPSO','_',func2str(f),'_Config'))
                        