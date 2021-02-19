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
addpath('Problems')

%% Configuration of data showing
fprint = 1; % Print data information each generations

% Problem parameters
NEI=1;
params.of = 2;% Objective functions
params.d=20; % Design variables 
params.c=0; % Constraints
params.np = 100; % Population size
params.gmax = 1000; % Max generations
params.lb=-15*ones(1,params.d); % Lower bounds
params.ub=30*ones(1,params.d); % Upper bounds

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
params.f = @Ackley;

%% Ejecucion de los algoritmos
parfor run = 1:NEI
    fprintf('************************* Corrida %d - GA *************************\n', run);
    GA(params, fprint, run);
end

%% Guardado de parametros
save(strcat('GA','_',func2str(params.f),'_Config'))