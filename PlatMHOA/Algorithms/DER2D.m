function DER2D(f, NP, G, D, DP, Fmin, Fmax, NEJ, xmin, xmax, fprint, run)
    
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
            
            r = randperm(NP,4);
            while(AP(r(1),D+1,k) > AP(r(2),D+1,k))
                r(1:2) = randperm(NP,2);
            end
            r(3:4) = randperm(NP,2);
            while(AP(r(3),D+1,k) > AP(r(4),D+1,k))
                r(3:4) = randperm(NP,2);
            end
            w = [AP(r(1),:,k); AP(r(2),:,k); AP(r(3),:,k); AP(r(4),:,k)];
   
            for j = 1:D
                APu(i,j) = w(1,j) + 0.5*F*(w(1,j) - w(2,j) + w(3,j) - w(4,j));
                if(APu(i,j) < xmin(j) || APu(i,j) > xmax(j))
                    APu(i,j) = xmin(j) + rand*(xmax(j) - xmin(j));
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
        
        if(mod(k,fprint) == 0)
            [~,idxbst] = min(AP(:,D+1,k));
            fprintf('G = %4d \t J* = %0.5f \t RV = %d \t x* = [', k, AP(idxbst,D+1,k), AP(idxbst,D+2,k))
            for j = 1:D-1, fprintf('%0.20f, ', AP(idxbst,j,k)), end, fprintf('%0.20f]\n', AP(idxbst,D,k));
        end
    end
    
    [~,~]  = mkdir('Data/DER2D/');
    OFNAME = func2str(f);
    namefile = strcat('Data/DER2D/','DER2D_',OFNAME,'_D',num2str(D),'_',num2str(run),'.mat');
    save(namefile)
    
end
