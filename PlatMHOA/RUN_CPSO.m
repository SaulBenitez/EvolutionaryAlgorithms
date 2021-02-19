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
addpath('Problems')

%% Parametros de impresion resultados
fprint = 1;

%% Parametros del problema
NEI = 1;
D = 20;
DP = 4;
NP = 50;
G = 1000;
NEJ = NP*G;
xmin = -15*ones(1,D);
xmax = 30*ones(1,D);

%% Parametros del algoritmo
C1 = 2.05;
C2 = 2.05;
w = 0.729;
vmin = -0.1*xmax;
vmax = 0.1*xmax;      
p = 0.7;
rho = 0.1;

%% Funcion objetivo
f = @Ackley;

%% Ejecucion de algoritmo
parfor run = 1:NEI
    fprintf('************************* Corrida %d - CPSO *************************\n', run);
    CPSO(f, NP, D, DP, G, NEJ, xmin, xmax, C1, C2, w, vmin, vmax, p, rho, fprint, run);
end

%% Guardado de parametros
save(strcat('CPSO','_',func2str(f),'_Config'))
                        