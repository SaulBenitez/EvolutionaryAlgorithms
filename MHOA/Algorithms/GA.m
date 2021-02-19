function GA(params, fprint, run)
    clc;
    
    AP = zeros(params.np, params.mofidx(end), params.gmax);
    
    %% Initialize population
    AP(:,:,1) = init_pop(params);
    
    %% Evaluate population
    for i=1:params.np
        AP(i,:,1)=evaluate_individual(AP(i,:,1),params);
    end
    
    [~,idxbst] = min(AP(:,params.d+1,1));
    fprintf('G = %4d \t J* = %0.5f \t RV = %d \t Je* = %0.15f \t Jev* = %5d \t x* = [', 1, AP(idxbst,params.d+1,1), AP(idxbst,params.d+2,1), AP(idxbst,params.d+3,1), AP(idxbst,params.d+4,1))
    for j = 1:params.d-1, fprintf('%0.20f, ', AP(idxbst,j,1)), end, fprintf('%0.20f]\n', AP(idxbst,params.d,1));
    
    %% Evolutionary cycle
    for gen=1:params.gmax-1
        %pop
        %% Crossover
        popu=crossover(AP(:,:,gen),params);
        %% Selection
        AP(:,:,gen+1) =selection(popu,params);
        
        if(mod(gen+1,fprint) == 0)
            [~,idxbst] = min(AP(:,params.d+1,gen+1));
            fprintf('G = %4d \t J* = %0.5f \t RV = %d \t Je* = %0.15f \t Jev* = %5d \t x* = [', gen+1, AP(idxbst,params.d+1,gen+1), AP(idxbst,params.d+2,gen+1), AP(idxbst,params.d+3,gen+1), AP(idxbst,params.d+4,gen+1))
            for j = 1:params.d-1, fprintf('%0.20f, ', AP(idxbst,j,gen+1)), end, fprintf('%0.20f]\n', AP(idxbst,params.d,gen+1));
        end
    end

    %% Print solution
    [~,~]  = mkdir('Data/GA/');
    OFNAME = func2str(params.f);
    namefile = strcat('Data/GA/','GA_',OFNAME,'_D',num2str(params.d),'_',num2str(run),'.mat');
    save(namefile)
end

%% Initialize population randomly in the search space
function pop=init_pop(params)
    pop = zeros(params.np, params.mofidx(end));

    for i = 1: params.np
        for j = 1:params.d
            pop(i,j) = params.lb(j) + rand()*(params.ub(j) - params.lb(j));
        end
    end
end

%% Evaluate an individual in population
function ind=evaluate_individual(ind,params)
    ind(params.ofidx:params.mofidx(end))=params.f(ind);
end

%% SBX crossover
function pop=crossover(pop,params)
    % Select parents by binary tournament
    candidate_indexes_a = [];
    candidate_indexes_b = [];
    parent_indexes = [];

    % Fill candidate vectors
    for i = 1:params.np
        candidate_indexes_a=[candidate_indexes_a;i];
        candidate_indexes_b=[candidate_indexes_b;i];
    end

    % Generate a random permutation
    candidate_indexes_a=shuffle(candidate_indexes_a);
    candidate_indexes_b=shuffle(candidate_indexes_b);

    % Binary comparision of candidates to select parents
    for i =1:params.np
        % If candidates have different constraint violation
        if pop(candidate_indexes_a(i),params.cvdidx) ~= pop(candidate_indexes_b(i),params.cvdidx)
            %Minimum constraint violation is selected
            if pop(candidate_indexes_a(i),params.cvdidx) < pop(candidate_indexes_b(i),params.cvdidx)
                parent_indexes=[parent_indexes;candidate_indexes_a(i)];
            else
                parent_indexes=[parent_indexes;candidate_indexes_b(i)];
            end
        else  % If candidates have the same constraint violation
            % Minimum objective function is selected
            if pop(candidate_indexes_a(i),params.ofidx) < pop(candidate_indexes_b(i),params.ofidx)
                parent_indexes=[parent_indexes;candidate_indexes_a(i)];
            else
                parent_indexes=[parent_indexes;candidate_indexes_b(i)];
            end
        end
    end

    % Generate offpsrings (SBX crossover and polynomial mutation)
    for i = 1:params.np/2
        % Use crossover probability
        if rand() < params.pc
            % Get the i-th couple
            parent_index_1 = parent_indexes(2 * (i-1) + 1);
            parent_index_2 = parent_indexes(2 * (i-1) + 2);

            offsprings = sbx(pop(parent_index_1,:), pop(parent_index_2,:), params);
            child_1=offsprings(1,:);
            child_2=offsprings(2,:);

            % Mutate offsprings
            child_1=poly_mutation(child_1,params);
            child_2=poly_mutation(child_2,params);

            % Evaluate offsprings
            child_1=evaluate_individual(child_1,params);
            child_2=evaluate_individual(child_2,params);

            % Add offsprings to the population
            pop=[pop;child_1];
            pop=[pop;child_2];
        end
    end
end

%% Random permutation
function v=shuffle(v)
    v=v(randperm(length(v)));
end

%% Replacement
function pop=selection(pop,params)
    % Sort regarding constraint violation and objective function value
    pop=sortrows(pop,[params.cvdidx params.ofidx]);
    % Remove worst solutions
    pop = pop(1:params.np,:);
end

%% SBX
function offspring=sbx(parent1, parent2, params)
    EPS = 1.0e-14;
    offspring=[parent1;parent2];

    for i = 1:params.d
        valueX1 = parent1(i);
        valueX2 = parent2(i);
        if rand() <= 0.5
            if (abs(valueX1 - valueX2) > EPS)
                if (valueX1 < valueX2)
                    y1 = valueX1;
                    y2 = valueX2;
                else
                    y1 = valueX2;
                    y2 = valueX1;
                end

                lowerBound = params.lb(i);
                upperBound = params.ub(i);

                ran = rand();
                beta = 1.0 + (2.0 * (y1 - lowerBound) / (y2 - y1));
                alpha = 2.0 - (beta)^(-(params.etac + 1.0));

                if (ran <= (1.0 / alpha))
                    betaq = (ran * alpha)^((1.0 / (params.etac + 1.0)));
                else
                    betaq = (1.0 / (2.0 - ran * alpha))^(1.0 / (params.etac + 1.0));
                end
                c1 = 0.5 * (y1 + y2 - betaq * (y2 - y1));

                beta = 1.0 + (2.0 * (upperBound - y2) / (y2 - y1));
                alpha = 2.0 - (beta)^(-(params.etac + 1.0));

                if (ran <= (1.0 / alpha))
                    betaq = ((ran * alpha))^((1.0 / (params.etac + 1.0)));
                else
                    betaq = (1.0 / (2.0 - ran * alpha))^( 1.0 / (params.etac + 1.0));
                end
                c2 = 0.5 * (y1 + y2 + betaq * (y2 - y1));

                if (c1 > params.ub(i))
                    c1 = params.ub(i);
                end
                if (c1 < params.lb(i))
                    c1 = params.lb(i);
                end

                if (c2 > params.ub(i))
                    c2 = params.ub(i);
                end
                if (c2 < params.lb(i))
                    c2 = params.lb(i);
                end


                if (rand() <= 0.5)
                    offspring(1,i) = c2;
                    offspring(2,i) = c1;
                else
                    offspring(1,i) = c1;
                    offspring(2,i) = c2;
                end
            else
                offspring(1,i) = valueX1;
                offspring(2,i) = valueX2;
            end
        else
            offspring(1,i) = valueX1;
            offspring(2,i) = valueX2;
        end
    end

end

%% Polynomial mutation
function ind=poly_mutation(ind,params)

    for i = 1:params.d
        if rand() <= params.pm
            y = ind(i);
            yl = params.lb(i);
            yu = params.ub(i);
            if yl == yu
                y = yl;
            else
                delta1 = (y - yl) / (yu - yl);
                delta2 = (yu - y) / (yu - yl);
                rnd = rand();
                mutPow = 1.0 / (params.etam + 1.0);
                if (rnd <= 0.5)
                    xy = 1.0 - delta1;
                    val = 2.0 * rnd + (1.0 - 2.0 * rnd) * ((xy)^( params.etam + 1.0));
                    deltaq = (val) ^( mutPow) - 1.0;
                else
                    xy = 1.0 - delta2;
                    val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * ((xy)^( params.etam + 1.0));
                    deltaq = 1.0 - (val)^(mutPow);
                end
                y = y + deltaq * (yu - yl);

                if y > params.ub(i)
                    y = params.ub(i);
                end
                if y < params.lb(i)
                    y = params.lb(i);
                end
            end
            ind(i) = y;
        end
    end
end

