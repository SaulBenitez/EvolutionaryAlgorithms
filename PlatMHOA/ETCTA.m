%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funcion objetivo: J = Je + Jev
% Sujeto a: xdot, Mp, P>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J] = ETCTA(p)

    %% Parametros de ponderacion 
    xi = [1000, 0.01];

    %% Parametros de tiempo
    to = 0; tf = 8; dt = 1e-3;
    ni = (tf-to)/dt + 1;
    t = linspace(to,tf,ni);
    tp = 2;
    nanf = 0;
    imagf = 1;

    %% Parametros dimencionales del sistema
    m = 3;
    n = 6;

    %% Parametros del control disparado por eventos
    rho = p(1:6);
    eps = p(7:9);
    Delta = diag(p(10:12));
    sigma = p(13);
    omega = p(14);
    psi = diag(p(15:17));
    
    u = zeros(m,ni);
    eventf = zeros(ni,1);
    
    %% Calculo de la matriz P
    P1 = zeros(m);
    P2 = zeros(m);
    P3 = zeros(m);
    for i = 1:m
        for j = 1:m
            if i==j
                P1(i,j) = sqrt( rho(i)*rho(i+m) + 2*sqrt( rho(i)^3/eps(i) ) );
                P2(i,j) = sqrt( rho(i+m)/eps(i) + sqrt( (4*rho(i))/eps(i)^3 ) );
                P3(i,j) = sqrt( rho(i)/eps(i) );
            end
        end
    end
    P = [P1, P3; P3, P2];
    
    %% Posiciones deseadas
    np = tf/tp;
    Qn = [pi/4, pi/4, pi/4, 0, 0, 0; pi, 0, 0, 0, 0, 0; pi/2, pi/4, pi/3, 0, 0, 0; 0, 0, 0, 0, 0, 0]';
    xd = zeros(n,ni);

    %% Inicializacion de variables de simulacion
    x = zeros(n,ni);
    z = zeros(n,ni);

    %% Ciclo principal
    if(isreal(p))
        cp = 1;
        uant = zeros(m,1);
        for i = 1:ni-1
            if(i*dt <= cp*tp+dt)
                xd(:,i) = Qn(:,cp);
            else 
                cp = cp + 1;
                xd(:,i) = Qn(:,cp);
            end

            [xdot, z(:,i), u(:,i), eventf(i)] = robot_rrr(x(:,i), xd(:,i), uant, P, Delta, sigma, omega, psi);
            uant = u(:,i);

            x(:,i+1) = x(:,i) + dt*xdot; 
            
            if (isnan(norm(x(:,i+1))) || ~isreal(x(:,i+1)))
                nanf = 1;
                break;
            end
        end
        imagf = 0;
    end
    
    if(nanf || imagf)
        J(1:3) = inf;
        J(4) = ni;
    else
        u(:,ni) = u(:,ni-1);
        xd(:,ni) = xd(:,ni-1);
        z(:,ni) = x(:,ni) - xd(:,ni);
        eventf(ni) = eventf(ni-1);
        
        %% Evaluacion de restricciones
        sizePP = (ni - 1)/np;
        countVR = 0;
%         if (det(P(1:4,1:4)) <= 0)
%            countVR = countVR + abs(det(P(1:4,1:4))); 
%         end
%         if (det(P(1:5,1:5)) <= 0)
%            countVR = countVR + abs(det(P(1:5,1:5))); 
%         end
%         if (det(P) <= 0)
%            countVR = countVR + abs(det(P)); 
%         end

        %% Evaluacion de funcion objetivo
        ITAE = zeros(1,m);
        for i = 1:m
            for j = 1:np
                if j == 1
                    ITAE(i) = ITAE(i) + sum(t(1:sizePP+1).*abs(z(i,1:sizePP+1)))*dt;
                else
                    ITAE(i) = ITAE(i) + sum(t(1:sizePP+1).*abs(z(i,sizePP*(j-1)+1:sizePP*j+1)))*dt;
                end
            end
        end
        
        J(4) = sum(eventf);
        J(3) = sum(ITAE);
        J(2) = countVR;
        J(1) = xi(1)*J(3) + xi(2)*J(4);

    end
    
end

function [xdot, z, u, eventf] = robot_rrr(x, xd, u, P, Delta, sigma, omega, psi)
    
    %% Matrices del sistema linelizado
    A = [zeros(3), eye(3); zeros(3,6)];
    B = [zeros(3); eye(3)];
    
    %% Inicializacion de matrices
    M = zeros(3);
    C = zeros(3);
    G = zeros(3,1);
    
    %% Parametros dinamicos del sistema
%     limu = 4.7124;
    limu = 5;
    l = [0.18, 0.15, 0.13];
    lc = [-0.0679573760623158050, 0.1186416627519466900, 0.0104239672085863310];
    Iz = [0.0718317628791155650, 0.0069993325747462372, 0.0013747544940748632];
    m = [5.8526474380885540000, 0.8769472564227329700, 0.3509389892379126700];
    Fv = diag([1.1697906594630387000, 0.6839555231774676600, 0.8017352022659838300]);

    g = 9.81; 
    
    M(1,1) = Iz(1) + Iz(2) + Iz(3) + l(1)^2*m(2) + l(1)^2*m(3) + l(2)^2*m(3) + lc(1)^2*m(1) + lc(2)^2*m(2) + lc(3)^2*m(3) + 2*l(1)*lc(3)*m(3)*cos(x(2) + x(3)) + 2*l(1)*l(2)*m(3)*cos(x(2)) + 2*l(1)*lc(2)*m(2)*cos(x(2)) + 2*l(2)*lc(3)*m(3)*cos(x(3));
    M(1,2) = m(3)*l(2)^2 + 2*m(3)*cos(x(3))*l(2)*lc(3) + l(1)*m(3)*cos(x(2))*l(2) + m(2)*lc(2)^2 + l(1)*m(2)*cos(x(2))*lc(2) + m(3)*lc(3)^2 + l(1)*m(3)*cos(x(2) + x(3))*lc(3) + Iz(2) + Iz(3);
    M(1,3) = Iz(3) + lc(3)^2*m(3) + l(1)*lc(3)*m(3)*cos(x(2) + x(3)) + l(2)*lc(3)*m(3)*cos(x(3));
    M(2,1) = M(1,2);
    M(2,2) = m(3)*l(2)^2 + 2*m(3)*cos(x(3))*l(2)*lc(3) + m(2)*lc(2)^2 + m(3)*lc(3)^2 + Iz(2) + Iz(3);
    M(2,3) = m(3)*lc(3)^2 + l(2)*m(3)*cos(x(3))*lc(3) + Iz(3);
    M(3,1) = M(1,3);
    M(3,2) = M(2,3);
    M(3,3) = m(3)*lc(3)^2 + Iz(3);
    
    C(1,1) = - 0;
    C(1,2) = - (2*x(4) + x(5) + 2*x(6))*l(1)*lc(3)*m(3)*sin(x(2) + x(3)) - (2*x(4) + x(5))*l(1)*l(2)*m(3)*sin(x(2)) - (2*x(4) + x(5))*l(1)*lc(2)*m(2)*sin(x(2)) - 2*x(6)*l(2)*lc(3)*m(3)*sin(x(3));
    C(1,3) = - (2*x(4) + x(6))*l(1)*lc(3)*m(3)*sin(x(2) + x(3)) - (2*x(4) + x(6))*l(2)*lc(3)*m(3)*sin(x(3));
    C(2,1) = x(4)*l(1)*lc(3)*m(3)*sin(x(2) + x(3)) + x(4)*l(1)*l(2)*m(3)*sin(x(2)) + x(4)*l(1)*lc(2)*m(2)*sin(x(2)) - 2*x(6)*l(2)*lc(3)*m(3)*sin(x(3));
    C(2,2) = 0;
    C(2,3) = - 2*x(5)*l(2)*lc(3)*m(3)*sin(x(3)) - x(6)*l(2)*lc(3)*m(3)*sin(x(3));
    C(3,1) = x(4)*l(2)*lc(3)*m(3)*sin(x(3)) + x(4)*l(1)*lc(3)*m(3)*sin(x(2) + x(3)) + 2*x(5)*l(2)*lc(3)*m(3)*sin(x(3));
    C(3,2) = x(5)*l(2)*lc(3)*m(3)*sin(x(3));
    C(3,3) = 0;
    
    G(1,1) = g*l(2)*m(3)*sin(x(1) + x(2)) + g*lc(2)*m(2)*sin(x(1) + x(2)) + g*l(1)*m(2)*sin(x(1)) + g*l(1)*m(3)*sin(x(1)) + g*lc(1)*m(1)*sin(x(1)) + g*lc(3)*m(3)*sin(x(1) + x(2) + x(3));
    G(2,1) = g*l(2)*m(3)*sin(x(1) + x(2)) + g*lc(2)*m(2)*sin(x(1) + x(2)) + g*lc(3)*m(3)*sin(x(1) + x(2) + x(3));
    G(3,1) = g*lc(3)*m(3)*sin(x(1) + x(2) + x(3));
    
    F = Fv*x(4:6);
    
    %% Control del robot
    z = x - xd;
%     V = 0.5*z'*P*z;
%     fx = [x(4:6); M\( - C*x(4:6) - G)];
%     gx = [zeros(3); M^-1];
%     ax = z'*P*fx;
%     bx = z'*P*gx;
    ax = z'*P*A*z;
    bx = z'*P*B;
    thbar = omega^2*bx*Delta*bx' - 2*ax*omega;
    
    ebar = - ax - bx*u - sigma*sqrt(ax^2 + thbar*bx*Delta*bx');
    if(ebar <= 0)
        if(norm(bx) == 0)
            gamma = 0;
        else
            gamma = (ax + sqrt(ax^2 + thbar*bx*Delta*bx'))/(bx*Delta*bx');
        end
%         u = - Delta*bx'*gamma;
%         u = - psi*bx'*gamma;
        u = - psi*Delta*bx'*gamma;
        eventf = 1;
    else
        eventf = 0;
    end
    
%     gamma = (ax + sqrt(ax^2 + thbar*bx*Delta*bx'))/(bx*Delta*bx');
%     u = - psi*Delta*bx'*gamma;
%     eventf = 1;
    
    for i = 1:3
        if(u(i) > limu)
            u(i) = limu;
        elseif(u(i) < -limu)
            u(i) = -limu;
        end
    end
    
    xdot = [x(4:6); M\(u - C*x(4:6) - G - F)];
    
end
    
