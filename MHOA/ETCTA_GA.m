%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autor:    Saul Benitez
% Titulo:   Genetic Algorithm (GA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------
% Description: This GA employs Binary Tournament (BT),
% Simulated Binary Crossover (SBX), Polynimial Mutation (PM)
% and Extintion Replacement.
% --------------------------------------------------------------   

clear; close all; clc;

% rng('shuffle');

%% Additional paths
addpath('Algorithms')

%% Configuration of data showing
fprint = 1; % Print data information each generations

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

% Problem parameters
NEI = 30;
params.of = 2;% Objective functions
params.d=17; % Design variables 
params.c=0; % Constraints
params.np = 50; % Population size
params.gmax = 2000; % Max generations
params.lb=[varrhoMin, epsMin, DeltaMin, sigmaMin, omegaMin, psiMin]; % Lower bounds
params.ub=[varrhoMax, epsMax, DeltaMax, sigmaMax, omegaMax, psiMax]; % Upper bounds

%% GA parameters     
params.pc = 0.9; % Crossover probability rate
params.pm = 0.6; % Mutation probability rate
params.etac = 20; % Crossover distribution index 
params.etam = 20; % Mutation distribution index 

% Matrix indexes
params.ofidx = params.d + 1; % Mono-bjective function index
params.cvdidx = params.ofidx + 1; % Constrain violation distance sum
params.mofidx = params.cvdidx+1:params.cvdidx+params.of;

% Objective function
params.f = @ETCTA;

%% Ejecucion de los algoritmos
parfor run = 1:NEI
    fprintf('************************* Corrida %d - GA *************************\n', run);
    GA(params, fprint, run);
end

%% Guardado de parametros
save(strcat('GA','_',func2str(f),'_Config'))