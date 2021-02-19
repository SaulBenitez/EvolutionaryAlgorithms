%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autor:    Saul Enrique Benitez Garcia                     %
% Archivo:  ParametricStatistics                            %
% Titulo:   Parametric Statistics of DE variants            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%********************************************************************
% Description: Parametric Statistics of DE Variants considers following 
% indexes:
%   - Mean(J): This is the mean of the mean of the population for each 
%   independent execution per algorithm.
%   - Standard Deviation(J): This is standard deviation of the mean for each 
%   independent execution per algorithm.
%   - Best(J): Best solution for every algorithm
%   - Worst(j): Worst solution for every algorith
%********************************************************************

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

%% Print in comand window parametric results based on the mean
[~, Meanidx] = min(POPavg);
[~, Stdidx] = min(POPstd);
[~, Bestidx] = min(BestMean);
[~, Worstidx] = min(WorstMean);
[~, BestAllidx] = min(BestAll(:,D+1));
[~, WorstAllidx] = min(BestAll(:,D+1));
for i = 1:NAE
    fprintf('\\textit{%s} & ', alg_name_tab{i})
    
    if(i == Meanidx)
        fprintf(' $\\mathbf{%0.4f}$ & ', POPavg(i))
    else
        fprintf(' $%0.4f$ & ', POPavg(i))
    end
    
    if(i == Stdidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  POPstd(i))
    else
        fprintf(' $%0.4f$ & ',  POPstd(i))
    end
    
    if(i == Bestidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  BestMean(i))
    else
        fprintf(' $%0.4f$ & ',  BestMean(i))
    end
    
    if(i == Worstidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  WorstMean(i))
    else
        fprintf(' $%0.4f$ & ',  WorstMean(i))
    end
    
    if(i == BestAllidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  BestAll(i,D+1))
    else
        fprintf(' $%0.4f$ & ',  BestAll(i,D+1))
    end
    
    if(i == WorstAllidx)
        fprintf(' $\\mathbf{%0.4f}$  \\\\\n',  WorstAll(i,D+1))
    else
        fprintf(' $%0.4f$ \\\\\n', WorstAll(i,D+1))
    end
    
end
fprintf('\n');
for i = 1:NAE
    fprintf('\\textit{%s} & ', alg_name_tab{i})
    
    if(i == Bestidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  BestMean(i))
    else
        fprintf(' $%0.4f$ & ',  BestMean(i))
    end
    
    if(i == Worstidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  WorstMean(i))
    else
        fprintf(' $%0.4f$ & ',  WorstMean(i))
    end
    
    if(i == BestAllidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  BestAll(i,D+1))
    else
        fprintf(' $%0.4f$ & ',  BestAll(i,D+1))
    end
    
    if(i == WorstAllidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  WorstAll(i,D+1))
    else
        fprintf(' $%0.4f$ & ', WorstAll(i,D+1))
    end
    
%     if(i == Meanidx)
%         fprintf(' $\\mathbf{%0.4f}$ & ', POPavg(i))
%     else
%         fprintf(' $%0.4f$ & ', POPavg(i))
%     end
    
    if(i == Stdidx)
        fprintf(' $\\mathbf{%0.4f}$ \\\\\n',  POPstd(i))
    else
        fprintf(' $%0.4f$ \\\\\n',  POPstd(i))
    end
    
end
fprintf('\n');

%%Impresion de resultados con base en solo el mejor
BestStd = zeros(NAE,1);
WorstBest = zeros(NAE,1);
BestSort = zeros(NEI,NAE);
for i = 1:NAE
    BestSort(:,i) = sort(POPbest(:,D+1,i), D+1);
    WorstBest(i) = BestSort(NEI,i);
    BestStd(i) = std(POPbest(:,D+1,i)); 
end
[~, WorstBestidx] = min(WorstBest);
[~, BestStdidx] = min(BestStd);
fprintf('\n');
for i = 1:NAE
    fprintf('\\textit{%s} & ', alg_name_tab{i})
    
    if(i == BestStdidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  BestStd(i))
    else
        fprintf(' $%0.4f$ & ',  BestStd(i))
    end
    
    if(i == Bestidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  BestMean(i))
    else
        fprintf(' $%0.4f$ & ',  BestMean(i))
    end
    
    if(i == BestAllidx)
        fprintf(' $\\mathbf{%0.4f}$ & ',  BestAll(i,D+1))
    else
        fprintf(' $%0.4f$ & ',  BestAll(i,D+1))
    end
    
    if(i == WorstBestidx)
        fprintf(' $\\mathbf{%0.4f}$ \\\\\n',  WorstBest(i))
    else
        fprintf(' $%0.4f$ \\\\\n', WorstBest(i))
    end
    
%     if(i == Stdidx)
%         fprintf(' $\\mathbf{%0.4f}$ \\\\\n',  POPstd(i))
%     else
%         fprintf(' $%0.4f$ \\\\\n',  POPstd(i))
%     end
    
end

%% Save optimal vector for each algorithm  
file = fopen('Best_Results.txt','w');
for i = 1:NAE
    [~,b] = min(POPbest(:,D+1,i));
    fprintf(file,'%s   \t\t      \n[',alg_name_tab{i});
    for j = 1:D-1
        fprintf(file,'%0.10f, ',POPbest(b,j,i));
    end
    fprintf(file,'%0.10f]\n',POPbest(b,D,i));
    fprintf(file,'J = %0.10f, R = %0.10f, J1 = %0.10f, J2 = %0.10f\n\n',POPbest(b,D+1,i), POPbest(b,D+2,i), POPbest(b,D+3,i), POPbest(b,D+4,i));
end
fclose(file);

%% Order and print mean
fprintf('============= Mean population Ordered =============\n')
[~, idx_sort] = sort(POPavg);
for i = 1:NAE
    fprintf('%s \t %f \n', alg_name_str{idx_sort(i)}, POPavg(idx_sort(i)))
end

%% Order and print STD
fprintf('============= Standard Deviation Ordered =============\n')
[~, idx_sort] = sort(POPstd);
for i = 1:NAE
    fprintf('%s \t %f \n', alg_name_str{idx_sort(i)}, POPstd(idx_sort(i)))
end

%% Order and print mean(Best) 
fprintf('============= mean(Best) Ordered =============\n')
[~, idx_sort] = sort(BestMean);
for i = 1:NAE
    fprintf('%s \t %f \n', alg_name_str{idx_sort(i)}, BestMean(idx_sort(i)))
end

%% Order and print Best
fprintf('============= Best Ordered =============\n')
[~, idx_sort] = sort(BestAll(:,D+1));
for i = 1:NAE
    fprintf('%s \t %f \n', alg_name_str{idx_sort(i)}, BestAll(idx_sort(i),D+1))
end

%% Order and print Worst of the Best
fprintf('============= Worst of the Best Ordered =============\n')
[~, idx_sort] = sort(WorstBest);
for i = 1:NAE
    fprintf('%s \t %f \n', alg_name_str{idx_sort(i)}, WorstBest(idx_sort(i)))
end

%% Order and print standard deviation for best samples
fprintf('============= Best standard deviation =============\n')
[~, idx_sort] = sort(BestStd);
for i = 1:NAE
    fprintf('%s \t %f \n', alg_name_str{idx_sort(i)}, BestStd(idx_sort(i)))
end

% %% Save data best with name for statistics
% for i = 1:NAE
%     file = fopen(alg_name_fil{i},'w');
%     for j = 1:NC
%         fprintf(file,'%0.20f\n', POPbest(j,D+1,i));
%     end
%     fclose(file);
% end

%% Graphics of convergence of algorithms
xwl = 10; ywl = 50; xwr = 1500; ywu = 700;
poswin = [xwl, ywl, xwr, ywu];
fsulab = 16;
fsuleg = 16;
fsaxesn = 12;
figure(1)
    set(gcf, 'Position', poswin)
    set(gca,'FontSize', fsaxesn)
    hold on
    plot(1:1:G, log10(POPhb(:,D+1,1)), '-.b');
    plot(1:1:G, log10(POPhb(:,D+1,2)), '--k');
    plot(1:1:G, log10(POPhb(:,D+1,3)), '-c');
    plot(1:1:G, log10(POPhb(:,D+1,4)), ':m');
    plot(1:1:G, log10(POPhb(:,D+1,5)), '-b');
    plot(1:1:G, log10(POPhb(:,D+1,6)), '-.r');
    plot(1:1:G, log10(POPhb(:,D+1,7)), '--b');
    plot(1:1:G, log10(POPhb(:,D+1,8)), ':g');
    plot(1:1:G, log10(POPhb(:,D+1,9)), '-k');
    plot(1:1:G, log10(POPhb(:,D+1,10)), '-.m');
    plot(1:1:G, log10(POPhb(:,D+1,11)), '-r');
    plot(1:1:G, log10(POPhb(:,D+1,12)), ':k');
    plot(1:1:G, log10(POPhb(:,D+1,13)), '--b');
    plot(1:1:G, log10(POPhb(:,D+1,14)), 'g');
    xlabel('$G$', 'Interpreter', 'latex', 'FontSize', fsulab, 'FontWeight', 'bold', 'Color', 'k')
    ylabel('$log(J)$', 'Interpreter', 'latex', 'FontSize', fsulab, 'FontWeight', 'bold', 'Color', 'k')
    legend({'\textit{DE/Rand/1/Bin}', '\textit{DE/Rand/1/Exp}', '\textit{DE/Best/1/Bin}', '\textit{DE/Best/1/Exp}', '\textit{DE/Current-to-Rand/1}', '\textit{DE/Current-to-Rand/1/Bin}', '\textit{DE/Current-to-Rand/1/Exp}', '\textit{DE/Current-to-Best/1}', '\textit{DE/Current-to-Best/1/Bin}', '\textit{DE/Current-to-Best/1/Exp}', '\textit{DE/Rand/2/Dir}', '\textit{DE/Rand/1/Either-or}', '\textit{GA}', '\textit{PSO}'}, 'Interpreter', 'latex', 'FontSize', fsuleg, 'Location', 'northeast')    
    % create a new pair of axes inside current figure
    axes('position',[0.2 0.5 .45 .25])
    box on % put box around new pair of axes
    hold on
    nsubx = 500;
    nsuby = 500;
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,1)), '-.b');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,2)), '--k');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,3)), '-r');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,4)), ':m');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,5)), '-b');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,6)), '-.r');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,7)), '--b');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,8)), ':g');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,9)), '-k');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,10)), '-.m');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,11)), '-r');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,12)), ':k');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,13)), '--b');
    plot(1:1:nsubx, log10(POPhb(1:nsuby,D+1,14)), 'g');
    axis([0,500,2.8,3.1])
    print -depsc convergencia.eps
