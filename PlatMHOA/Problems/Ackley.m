function J = Ackley(x)
    
    % Ackley Function -> f(0, 0,...,0) = 0
    n = 20;
    a = 20; b = 0.2; c = 2*pi;
    s1 = 0; s2 = 0;
    for i=1:n
       s1 = s1+x(i)^2;
       s2 = s2+cos(c*x(i));
    end
    f = -a*exp(-b*sqrt(1/n*s1))-exp(1/n*s2)+a+exp(1);
    g = 0;
    
    J(1) = f;
    J(2) = g;
    J(3) = 0;
    J(4) = 0;
end
