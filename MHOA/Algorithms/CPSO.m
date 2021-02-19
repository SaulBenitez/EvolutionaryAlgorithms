%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autor:    Saul Enrique Benitez Garcia                     %
% Creado:   5/Febrero/2020                                   %
% Archivo:  pso.m                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CPSO(f, NP, D, DP, G, NEJ, xmin, xmax, C1, C2, w, vmin, vmax, p, rho, fprint, run)
    
    Q = 300; % Exploitation iterations
    countJE = 0; %Objective function evaluations counter
    
    X = zeros(NP,D+DP,G); % Particles position 
    V = zeros(NP,D,G); % Particles Speed
    Lbest = zeros(NP,D+DP,G); % Local Best postion per particle 
    Gbest = zeros(G,D+DP); % Global Best swarm position 
    
    k = 1; % Generation counter
    
    %/////////////////////////////////////////////////////////
    % Initialization first particle swarm (1st generation)
    %/////////////////////////////////////////////////////////
    X(:,1:D,k) = xmin + rand(NP,D).*(xmax - xmin);
    V(:,1:D,k) = vmin + rand(NP,D).*(vmax - vmin);
    for i = 1:NP
        X(i,D+1:D+DP,k) = f(X(i,1:D,k));
    end
    countJE = countJE + NP;
    
    %///////////////////////////////////////////////////
    % Get best particle swarm and best particle history
    %///////////////////////////////////////////////////
    Lbest(:,:,k) = X(:,:,k);
    [~,gbstidx] = min(X(:,D+1,k));
    Gbest(k,:) = Lbest(gbstidx,:,k);
    
    
    %//////////////////////////
    % Print Results
    %//////////////////////////
%     fprintf('G = %4d \t J* = %0.5f \t RV = %d \t Je* = %0.15f \t Jev* = %5d \t x* = [', k, Gbest(k,D+1), Gbest(k,D+2), Gbest(k,D+3), Gbest(k,D+4))
    fprintf('G = %4d \t J* = %0.5f \t RV = %d \t x* = [', k, Gbest(k,D+1), Gbest(k,D+2))
    for j = 1:D-1, fprintf('%0.20f, ', Gbest(k,j)), end, fprintf('%0.20f]\n',Gbest(k,D));
    
    
    %//////////////////////////////////
    % PSO main loop
    %//////////////////////////////////
    r = (xmax - xmin)/2; % Set initial chaotic search radious
    while(k < G && countJE < NEJ)
%         Favg = mean(X(:,D+1,k));
        for i = 1:NP 
%             %+++++++++++++++++++++++++++++++++++++++++++ 
%             % Adaptive inertia weight factor
%             %+++++++++++++++++++++++++++++++++++++++++++ 
%             if(X(i,D+1,k) <= Favg)
%                 w = wmin + eta*(X(i,D+1,k) - Lbest(i,D+1,k))/(Favg - Lbest(i,D+1,k)); 
%             else
%                 w = wmax;
%             end
%             %+++++++++++++++++++++++++++++++++++++++++++ 
%             % Linear decrease weight factor
%             %+++++++++++++++++++++++++++++++++++++++++++
%             w = wmax - (wmax - wmin)*k/G;
%             
%             %+++++++++++++++++++++++++++++++++++++++++++
%             % Time-Varying Inertia Weight
%             %+++++++++++++++++++++++++++++++++++++++++++
%             w = (wo - wf)*(G - k)/k + wf;
            
            %+++++++++++++++++++++++++++++++++++++++++++ 
            % Update speed and position particles
            %+++++++++++++++++++++++++++++++++++++++++++ 
%             V(i,:,k+1) = w*V(i,:,k) + C1*rand()*(Lbest(i,1:D,k) -  X(i,1:D,k)) + C2*rand()*(Gbest(k,1:D) -  X(i,1:D,k));
            V(i,:,k+1) = w*(V(i,:,k) + C1*rand()*(Lbest(i,1:D,k) -  X(i,1:D,k)) + C2*rand()*(Gbest(k,1:D) -  X(i,1:D,k)));
            X(i,1:D,k+1) = X(i,1:D,k) + V(i,:,k+1);
            
            %+++++++++++++++++++++++++++++++++++++++++++    
            % Apply boundary constraint-handling methods
            %+++++++++++++++++++++++++++++++++++++++++++
            for j = 1:D
                if(X(i,j,k+1) < xmin(j))
%                     PS(i,j,k+1) = xmin(j); % Boundary method
                    X(i,j,k+1) = xmin(j) + rand()*(xmax(j) - xmin(j)); % Random method
%                     PS(i,j,k+1) = 2*xmin(j) - PS(i,j,k+1); % Reflection
%                     PS(i,j,k+1) = xmax(j) - mod(xmin(j) - PS(i,j,k+1),abs(xmax(j) - xmin(j))); % Wrapping
%                     alfa = rand();
%                     PS(i,j,k+1) = alfa*xmin(j) + (1 - alfa)*PS(i,j,k); % Evolutionary
%                     V(i,j,k+1) = 0; 
                elseif(X(i,j,k+1) > xmax(j))
%                     PS(i,j,k+1) = xmax(j); % Boundary method
                    X(i,j,k+1) = xmin(j) + rand()*(xmax(j) - xmin(j)); % Random method
%                     PS(i,j,k+1) = 2*xmax(j) - PS(i,j,k+1); % Reflection
%                     PS(i,j,k+1) = xmin(j) - mod(PS(i,j,k+1) - xmax(j),abs(xmax(j) - xmin(j))); % Wrapping
%                     beta = rand();
%                     PS(i,j,k+1) = beta*xmax(j) + (1 - beta)*Pbest(i,j,k); % Evolutionary
%                     V(i,j,k+1) = 0; 
                end
            end

            %+++++++++++++++++++++++++++++++++++++++
            % Evaluate the new position performance
            %+++++++++++++++++++++++++++++++++++++++     
            X(i,D+1:D+DP,k+1) = f(X(i,1:D,k+1));
            
            %+++++++++++++++++++++++++++++++++++++++
            % Update local best position
            %+++++++++++++++++++++++++++++++++++++++  
            if(X(i,D+1,k+1) < Lbest(i,D+1,k))
                Lbest(i,:,k+1) = X(i,:,k+1);
            else
                Lbest(i,:,k+1) = Lbest(i,:,k);
            end
        end
        
        countJE = countJE + NP;
        k = k + 1;
        
        %///////////////////////////
        % Get best particle swarm 
        %///////////////////////////
        [~,gbstidx] = min(Lbest(:,D+1,k));
        Gbest(k,:) = Lbest(gbstidx,:,k);
                                                                
        %///////////////////////////////////////////////////////
        % Apply Chaotic Local Search for Global best particle
        %///////////////////////////////////////////////////////
        z = Gbest(k,:); % Original global best
        y = z*1.1; % Exploitation results
        i = 1;
        flag = 1;
        while(flag && i <= Q)
            cx = rand(1,D);
            for j = 1:D
                if(cx(j) < p)
                    cx(j) = cx(j)/p;
                else
                    cx(j) = (1 - cx(j))/(1 - p);
                end
            end
            y(1:D) = z(1:D) + r.*(2*cx - 1);
            y(D+1:D+DP) = f(y(1:D));
            countJE = countJE + 1;
            if(y(D+1) < z(D+1))
                z = y;
                flag = 0;
            end
            i = i + 1;
        end
        Gbest(k,:) = z;  
        r = rho*r;
        
        %//////////////////////////
        % Print Results
        %//////////////////////////
        if(mod(k,fprint) == 0)
%             fprintf('G = %4d \t J* = %0.5f \t RV = %d \t Je* = %0.15f \t Jev* = %5d \t x* = [', k, Gbest(k,D+1), Gbest(k,D+2), Gbest(k,D+3), Gbest(k,D+4))
            fprintf('G = %4d \t J* = %0.5f \t RV = %d \t x* = [', k, Gbest(k,D+1), Gbest(k,D+2))
            for j = 1:D-1, fprintf('%0.20f, ', Gbest(k,j)), end, fprintf('%0.20f]\n',Gbest(k,D));
        end
        
    end
    AP = Lbest;
    for i = k+1:G
        AP(:,:,i) = Lbest(:,:,k);
        Gbest(i,:) = Gbest(k,:);
    end
    for i = 1:G
        AP(1,:,i) = Gbest(i,:);
    end
    
    [~,~]  = mkdir('Data/CPSO/');
    OFNAME = func2str(f);
    namefile = strcat('Data/CPSO/','CPSO_',OFNAME,'_D',num2str(D),'_',num2str(run),'.mat');
    save(namefile)
end


% function [x] = cls(x, D, DP, xmin, xmax) % Chaotic Local Search
%     
%     iota = xmax - xmin; 
%     y = x*1.1; % Exploitation results
%     explk = 300; % Exploitation iterations
%     
%     %% Chaotic local search main loop
%     k = 1;
%     flag = 1;
%     while(flag && k <= explk)
%         cx = (x(1:D) -  xmin)/iota;
%         cx = 4*cx.*(1 - cx);
%         y(1:D) = xmin + cx*iota;
%         y(D+1:D+DP) = f(y(1:D));
%         if(y(D+1) < x(D+1))
%             x = y;
%             flag = 1;
%         else
%             flag = 0;
%         end
%         k = k + 1;
%     end
%     
% end
