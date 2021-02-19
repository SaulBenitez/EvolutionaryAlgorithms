function J=problem(x)
    %% Beale function -> f(3,0.5)=0
    % f = (1.5-x(1)+x(1)*x(2))^2 + (2.25-x(1)+x(1)*x(2)^2)^2 + (2.625 -x(1) +x(1)*x(2)^3)^2;
    % g = 0;
    % res = [f ; g];

    %% Rosenbrock function -> f(1,1)=0
    % f = (1-x(1))^2 + 100*(x(2)-x(1)^2)^2;
    % g1 = (x(1)-1)^3 - x(2) + 1;
    % g2 = x(1) + x(2) -2;
    % res = [f ; g1 ; g2];

    % %% G01 function -> f(1,1,1,1,1,1,1,1,1,3,3,3,1)=-15
    % sum1=0; 
    % sum2=0; 
    % sum3=0;
    % for i=1:4
    %     sum1=sum1+x(i);
    % end
    % for i=1:4
    %     sum2=sum2+x(i)^2;
    % end
    % for i=5:13
    %     sum3=sum3+x(i);
    % end
    % f = 5*sum1 - 5*sum2 -sum3;
    % g1 = 2*x(1) + 2*x(2) + x(10) + x(11) - 10;
    % g2 = 2*x(1) + 2*x(3) + x(10) + x(11) - 10;
    % g3 = 2*x(2) + 2*x(3) + x(11) + x(12) - 10;
    % g4 = -8*x(1) + x(10);
    % g5 = -8*x(2) + x(11);
    % g6 = -8*x(3) + x(12);
    % g7 = -2*x(4) - x(5) + x(10);
    % g8 = -2*x(6) - x(7) + x(11);
    % g9 = -2*x(8) - x(9) + x(12);
    % res = [f; g1; g2; g3; g4; g5; g6; g7; g8; g9];

%     %% Beale Function -> f(3, 0.5) = 0
%     f = (1.5-x(1)*(1-x(2)))^2+(2.25-x(1)*(1-x(2)^2))^2+(2.625-x(1)*(1-x(2)^3))^2;
%     g = 0;
% 
%     J(1) = f;
%     J(2) = g;
%     J(3) = 0;
%     J(4) = 0;
    
    % Ackley Function -> f(0, 0,...,0) = 0
    n = 2;
    a = 20; b = 0.2; c = 2*pi;
    s1 = 0; s2 = 0;
    for i=1:n
       s1 = s1+x(i)^2;
       s2 = s2+cos(c*x(i));
    end
    f = -a*exp(-b*sqrt(1/n*s1))-exp(1/n*s2)+a+exp(1);
    g = 0;

%     % Bohachecsky function 1 -> f(0, 0,...,0) = 0
%     % The number of variables n = 2. -100<= xi <= 100
%     f = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))-0.4*cos(4*pi*x(2))+0.7;
%     g = 0;

%     % Branin function -> f(0, 0) = 0
%     % The number of variables n = 2. -5 <= x1 <= 10, 0 <= x2 <= 15
%     f = (x(2)-(5.1/(4*pi^2))*x(1)^2+5*x(1)/pi-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
%     g = 0;

%     % Colville function -> f(0, 0, 0, 0) = 0
%     % The number of variables n = 4.
%     f  = 100*(x(1)^2-x(2))^2+(x(1)-1)^2+(x(3)-1)^2+90*(x(3)^2-x(4))^2+...
%     10.1*((x(2)-1)^2+(x(4)-1)^2)+19.8*(x(2)^-1)*(x(4)-1);
%     g = 0;

%     % Easom function -> f(pi, pi) = 0
%     % The number of variables n = 2.
%     f = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
%     g = 0;

%     % Ackley function -> f(pi, pi) = 0
%     % The number of variables n = 2.
%     n = 200;
%     f = -20*exp(-0.2*sqrt(sum(x.^2)/n)) - exp(sum(cos(2*pi*x))/n) + 20 + exp(1);
%     g = 0;
    
    
    J(1) = f;
    J(2) = g;
    J(3) = 0;
    J(4) = 0;
end
