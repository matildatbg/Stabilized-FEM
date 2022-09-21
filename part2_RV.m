

%%% RV-method
clear all; close all;

x1_0 = 0.3;
x2_0 = 0;
r0 = 0.25;
CFL = 0.5;
C_vel = 0.25;
C_RV = 1;
fprime = @(x1,x2) 2*pi*[-x2; x1];
hmax = [1/4 1/8 1/16 1/32];
T = 1;

% Choose type of IC
n = input('\nType \n 1 for Problem 1.1 \n 2 For Problem 1.3:   ');
switch n
    case 1
        u0 = @(x1,x2) 1/2.*(1-tanh(((x1-x1_0).^2+(x2-x2_0).^2)./(r0^2)-1));
        
    case 2
        u0 = @(x1,x2) ((x1-x1_0).^2+(x2-x2_0).^2)<=r0^2 *1;
end

for i=1:length(hmax)
    
    % Setting domain
    geometry = @circleg;
    [p,e,t] = initmesh(geometry, 'hmax', hmax(i));
    x1 = p(1,:);
    x2 = p(2,:);
    numElem = size(t,2);
    
    u_initial = double(u0(x1,x2))';
    convTerm = 2*pi*[-x2; x1];
    
    % Stable time step
    kn = CFL*hmax(i)/max(sqrt(convTerm(1,:).^2 + convTerm(2,:).^2));
    kn = T/ceil(T/kn);
    
    % Assemble matrices
    M = mass2D(p,t);
    C = convMat2D(p, t, convTerm(1,:), convTerm(2,:));
    
    crankNic = @(Rv) (2*M+kn*C+kn*Rv) \ (2*M-kn*C-kn*Rv);
    I = speye(length(p));
    
    % Timestepping
    time = 0;
    xi = u_initial;
    xi_prev = xi;
    while time <= T
        
        % calculating and normalizing residual
        Ruh = M\(1/kn*( M*xi-M*xi_prev) + C*(xi));
        Ruh = Ruh/max((xi-mean(xi)));
        
        % calculating the artificial viscosity based on the residual
        for K = 1:numElem
            nodes = t(1:3,K);
            coords = p(:,nodes);
            d = [norm(coords(:,1)-coords(:,2)) norm(coords(:,2)-coords(:,3))...
                norm(coords(:,3)-coords(:,1))];
            h_K = min(d);
            
            fp = fprime(coords(1,:), coords(2,:));
            
            beta_K = max(sqrt(fp(1,:).^2+fp(2,:).^2));
            Res_K = max(abs(Ruh(nodes)));
            eps(K) = min(C_vel*h_K*beta_K, C_RV*h_K^2*Res_K);
        end
        
        % Assembling matrices
        Rv = RvMat2D(p,t,eps);
        A = sparse(crankNic(Rv));
        A(e(1 ,:) ,:) = I(e(1 ,:) ,:);
        
        xi_prev = xi;
        
        xi = A*xi_prev;
        xi(e(1,:)) = 0;
        % Draws every update
        %pdesurf(p,t,xi);
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
    a = pdesurf(p,t,xi);
    %%a.EdgeColor = [0 0 0];
    title(["u_h", "h_{max} = " + hmax(i) + ", T = " + T]);
    xlabel("x_1");
    ylabel("x_2");
    zlabel("u_h(x_1, x_2)");
    
    colormap('default');
    
    error = u_initial-xi;
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




