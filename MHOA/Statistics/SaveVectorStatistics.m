%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autor:    Saul Enrique Benitez Garcia                     %
% Archivo:  ParametricStatistics                            %
% Titulo:   Parametric Statistics of DE variants            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%================== Data samples in PSO algorithm =======================
clear; close all; format short; clc;
load('ETC-AlgorithmsDataSamples')

%% Names for files and tables
alg_name_fil = {'DE Rand 1 Bin'; 
            'DE Rand 1 Exp'; 
            'DE Best 1 Bin'; 
            'DE Best 1 Exp'; 
            'DE Current-to-Rand 1'; 
            'DE Current-to-Rand 1 Bin'; 
            'DE Current-to-Rand 1 Exp'; 
            'DE Current-to-Best 1'; 
            'DE Current-to-Best 1 Bin'; 
            'DE Current-to-Best 1 Exp'; 
            'DE Rand 2 Dir'; 
            'DE Rand 1 Either-or';
            'GA';
            'PSO';
            'DE-PSO'};

alg_name_tab = {'DE/Rand/1/Bin'; 
            'DE/Rand/1/Exp'; 
            'DE/Best/1/Bin'; 
            'DE/Best/1/Exp'; 
            'DE/Current-to-Rand/1'; 
            'DE/Current-to-Rand/1/Bin'; 
            'DE/Current-to-Rand/1/Exp'; 
            'DE/Current-to-Best/1'; 
            'DE/Current-to-Best/1/Bin'; 
            'DE/Current-to-Best/1/Exp'; 
            'DE/Rand/2/Dir'; 
            'DE/Rand/1/Either-or';
            'GA';
            'PSO';
            'DE-PSO'};
        
alg_name_str = ["DE/Rand/1/Bin              ",...
                "DE/Rand/1/Exp              ",... 
                "DE/Best/1/Bin              ",... 
                "DE/Best/1/Exp              ",... 
                "DE/Current-to-Rand/1       ",... 
                "DE/Current-to-Rand/1/Bin   ",... 
                "DE/Current-to-Rand/1/Exp   ",... 
                "DE/Current-to-Best/1       ",...
                "DE/Current-to-Best/1/Bin   ",...
                "DE/Current-to-Best/1/Exp   ",...
                "DE/Rand/2/Dir              ",...
                "DE/Rand/1/Either-or        ",...
                "GA                         ",...
                "PSO                        ",...
                "DE-PSO                     "];

%% Get best mean and worst            
POPsort = zeros(size(POP));
POPIDXsort = zeros(NP,NEI,NAE);
POPmean = zeros(NEI,NAE);
POPstd = zeros(NAE,1);
POPavg = zeros(NAE,1);
POPbest = zeros(NEI,D+DP,NAE);
POPworst = zeros(NEI,D+DP,NAE);
BestMean = zeros(NAE,1);
WorstMean = zeros(NAE,1);
BestAll = zeros(NAE,D+DP);
WorstAll = zeros(NAE,D+DP);
for i = 1:NAE
    for j = 1:NEI
        [POPsort(:,:,j,i), POPIDXsort(:,j,i)] = sortrows(POP(:,:,j,i), [D+2,D+1]);
        POPmean(j,i) = mean(POPsort(:,D+1,j,i));
        POPbest(j,:,i) = POPsort(1,:,j,i);
        POPworst(j,:,i) = POPsort(NP,:,j,i);
    end
    POPstd(i) = std(POPmean(:,i)); 
    POPavg(i) = mean(POPmean(:,i)); 
    BestMean(i) = mean(POPbest(:,D+1,i));
    WorstMean(i) = mean(POPworst(:,D+1,i));
    [~,b] = min(POPbest(:,D+1,i));
    BestAll(i,:) = POPbest(b,:,i);
    [~,b] = max(POPworst(:,D+1,i));
    WorstAll(i,:) = POPworst(b,:,i);
end


%% Save data best with name for statistics
for i = 1:NAE
    file = fopen(['./groups/', alg_name_fil{i}],'w');
    for j = 1:NEI
        fprintf(file,'%0.20f\n', POPbest(j,D+1,i));
    end
    fclose(file);
end
% fclose(file);
