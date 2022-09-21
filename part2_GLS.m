
%%% FEM Part 2 - Stabilization with Galerkin Least Squares

clear all; close all;

x1_0 = 0.3;
x2_0 = 0;
r0 = 0.25;
CFL = 0.5;
hmax = [1/4 1/8 1/16 1/32];
T = 1;
%delta = 0.1*hmax;

% Choose stuff
n = input('\nType \n 1 for smooth IC \n 2 for step IC:   ');
switch n
    case 1
        u0 = @(x1,x2) 1/2.*(1-tanh(((x1-x1_0).^2+(x2-x2_0).^2)./(r0^2)-1));
        
    case 2
        u0 = @(x1,x2) ((x1-x1_0).^2+(x2-x2_0).^2)<=r0^2 *1;
end

for i=1:length(hmax)
    geometry = @circleg;
    [p,e,t] = initmesh(geometry, 'hmax', hmax(i));
    x1 = p(1,:);
    x2 = p(2,:);
    u_initial = double(u0(x1,x2))';
    convTerm = 2*pi*[-x2; x1];
    delta = (1/2)*hmax(i)/max(sqrt(convTerm(1,:).^2 + convTerm(2,:).^2));
    
    
    % Assemble matrices
    M = mass2D(p,t);
    C = convMat2D(p, t, convTerm(1,:), convTerm(2,:));
    C_tran = C';
    Sd = SDMat2D(p, t, convTerm(1,:), convTerm(2,:));
    
    % Stable time step
    kn = CFL*hmax(i)/max(sqrt(convTerm(1,:).^2 + convTerm(2,:).^2));
    kn = T/ceil(T/kn);
    
    % Crank-Nicolson:
    nom = 1/kn*(M+delta*C_tran) - 1/2*(C+delta*Sd);
    denom = 1/kn*(M+delta*C_tran) + 1/2*(C+delta*Sd);
    
    A = denom\nom;
    
    % Sets boundary condition
    I = eye(length(p));
    A(e(1 ,:) ,:) = I(e(1 ,:) ,:);
    
    uh = u_initial;
    %pdesurf(p,t,uh);
    
    time = 0;
    while (time <= T)
        %set boundary strongly
        uh(e(1,:)) = 0;
        uh = A*uh;
        % Draws every update
        %pdesurf(p,t,uh);
        %drawnow;
        
        time = time+kn;
    end
    
    figure(i);
    subplot(2,1,1)
    pdesurf(p,t,u_initial);
    title(["Exact solution u", "h_{max} = " + hmax(i) + ", T = " + T]);
    xlabel("x_1");
    ylabel("x_2");
    zlabel("u(x_1, x_2)");
    
    subplot(2,1,2)
    a = pdesurf(p,t,uh);
    %%a.EdgeColor = [0 0 0];
    title(["u_h", "h_{max} = " + hmax(i) + ", T = " + T]);
    xlabel("x_1");
    ylabel("x_2");
    zlabel("u_h(x_1, x_2)");
    
    colormap('default');
    
    error = u_initial-uh;
    L2E(i) = sqrt(error'*M*error);
end

P = polyfit(log10(hmax), log10(L2E),1);
figure(i+1)
loglog(hmax,L2E,'-b', hmax, hmax.^(P(1)), '-r');
title("Error");
xlabel("log h_{max}")
ylabel("log ||u-u_h||_{L^2(\Omega)}");
legend("Error", "h_{max}^P, P =" + P(1),'Location', 'southeast');

%Convergence between each hmax
q = zeros(1,length(hmax)-1);
for i=1:length(q)
    q(i) = (log10(L2E(i+1))-log10(L2E(i)))/(log10(hmax(i+1))-log10(hmax(i)));
end

