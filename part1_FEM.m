
%%% Solves the linear advection equation 
%%% with GFEM using two different initial
%%% conditions: one smooth gaussian function
%%% and one circular step function.
%%%
%%% Calls the function mass2D.m and convMat2D.m

clear all; close all;

x1_0 = 0.3;
x2_0 = 0;
r0 = 0.25;
CFL = 0.5;
hmax = [1/4 1/8 1/16 1/32];
T = 1;

% Choose initial condition
n = input('\nType \n 1 for Problem 1.1 \n 2 For Problem 1.3:   ');
switch n
    case 1
        u0 = @(x1,x2) 1/2.*(1-tanh(((x1-x1_0).^2+(x2-x2_0).^2)./(r0^2)-1));
        directory = "./Part_1_1/";
        mkdir(directory);
    case 2
        u0 = @(x1,x2) ((x1-x1_0).^2+(x2-x2_0).^2)<=r0^2 *1;
        directory = "./Part_1_3/";
        mkdir(directory);
end

for i=1:length(hmax)
    % sets geometry
    geometry = @circleg; 
    [p,e,t] = initmesh(geometry, 'hmax', hmax(i));

    % sets values
    x1 = p(1,:);
    x2 = p(2,:);
    u_initial = double(u0(x1,x2))';
    convTerm = 2*pi*[-x2; x1];

    % Assemble matrices
    M = mass2D(p,t);
    C = convMat2D(p, t, convTerm(1,:), convTerm(2,:));
    
    % Stable time step
    kn = CFL*hmax(i)/max(sqrt(convTerm(1,:).^2 + convTerm(2,:).^2)); 
    kn = T/ceil(T/kn);
    
    % Crank Nicolson
    nom = 2*M-kn*C;
    denom = 2*M+kn*C;

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
    
    % Save the figures
    filename = directory +"plot_h_" + num2str(hmax(i)) + ".png"; 
    saveas(gcf, filename);
end


P = polyfit(log10(hmax), log10(L2E),1);
figure(i+1)
loglog(hmax,L2E,'-b', hmax, hmax.^(P(1)), '-r');
title("Error");
xlabel("log h_{max}")
ylabel("log ||u-u_h||_{L^2(\Omega)}");
legend("Error", "h_{max}^P, P =" + P(1),'Location', 'southeast');
saveas(gcf, directory+'Error.png');

%Convergence between each hmax
q = zeros(1,length(hmax)-1);
for i=1:length(q)
    q(i) = (log10(L2E(i+1))-log10(L2E(i)))/(log10(hmax(i+1))-log10(hmax(i)));
end


