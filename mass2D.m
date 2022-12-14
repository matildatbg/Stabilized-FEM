
%%% Copied from the book
%%% "The Finite Element Method: Theory, Implementation, and
%%% Practice%: Larson, Bengtzon 
%%%
%%% Creates mass matrix for 2D-problems.

function M = mass2D(p,t)
N = size(p,2);
M = sparse(N,N);
%M = zeros(N,N);
for K =1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    area_K = polyarea(x,y);
    
    Mk = (1/12)*[2 1 1; 1 2 1; 1 1 2]*area_K;
    M(nodes, nodes) = M(nodes,nodes) + Mk;
end
