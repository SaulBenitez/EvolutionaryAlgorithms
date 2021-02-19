function DECB1B(f, NP, G, D, DP, Fmin, Fmax, Kmin, Kmax, CR, NEJ, xmin, xmax, fprint, run)
    
    countJE = 0;
    AP = zeros(NP,D+DP,G);
    APu = zeros(NP,D+DP);

    k = 1; 
    
    for i = 1:NP
        AP(i,1:D,k) = xmin + rand(1,D).*(xmax - xmin);
        AP(i,D+1:D+DP,k) = f(AP(i,1:D,k));
    end
    
    countJE = countJE + NP;
    
    [~,idxbst] = min(AP(:,D+1,k));
    fprintf('G = %4d \t J* = %0.5f \t RV = %d \t x* = [', k, AP(idxbst,D+1,k), AP(idxbst,D+2,k))
    for j = 1:D-1, fprintf('%0.20f, ', AP(idxbst,j,k)), end, fprintf('%0.20f]\n', AP(idxbst,D,k));
  
    while(k < G && countJE <= NEJ)
        for i = 1:NP
            F = Fmin + rand()*(Fmax - Fmin);
            K = Kmin + rand()*(Kmax - Kmin);
            
            r = randperm(NP,4);
            r(i==r) = [];

            jrand = randi(D);
            for j = 1:D
                if(rand < CR || j == jrand)
                    APu(i,j) = AP(i,j,k) + K*(AP(idxbst,j,k) - AP(i,j,k)) + F*(AP(r(1),j,k) - AP(r(2),j,k));
                    if(APu(i,j) < xmin(j) || APu(i,j) > xmax(j))
                        APu(i,j) = xmin(j) + rand(1,1)*(xmax(j) - xmin(j));
                    end
                else
                    APu(i,j) = AP(i,j,k);
                end
            end
  
            APu(i,D+1:D+DP) = f(APu(i,1:D));
            if(APu(i,D+1) < AP(i,D+1,k))
                AP(i,:,k+1) = APu(i,:);
            else
                AP(i,:,k+1) = AP(i,:,k);
            end
 
        end
        
        countJE = countJE + NP;
        k = k + 1;
        
        [~,idxbst] = min(AP(:,D+1,k));
        if(mod(k,fprint) == 0)
            fprintf('G = %4d \t J* = %0.5f \t RV = %d \t x* = [', k, AP(idxbst,D+1,k), AP(idxbst,D+2,k))
            for j = 1:D-1, fprintf('%0.20f, ', AP(idxbst,j,k)), end, fprintf('%0.20f]\n', AP(idxbst,D,k));
        end
    end
    
    [~,~]  = mkdir('Data/DECB1B/');
    OFNAME = func2str(f);
    namefile = strcat('Data/DECB1B/','DECB1B_',OFNAME,'_D',num2str(D),'_',num2str(run),'.mat');
    save(namefile)
    
end
