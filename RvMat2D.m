%%% Copied from the book
%%% "The Finite Element Method: Theory, Implementation, and
%%% Practice%: Larson, Bengtzon
%%% with a slight modificaiton to introduce epsilon as stabilization.
%%%
%%% Constructs stabilization for 2D-problems.
function Rv = RvMat2D(p,t, eps)

N = size(p,2);
Rv = sparse(N,N);
%Rv = zeros(N,N);

for K =1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    area_K = polyarea(x,y);
    bi = [y(2)-y(3); y(3)-y(1); y(1)-y(2)]/(2*area_K);
    ci = [x(3)-x(2); x(1)-x(3); x(2)-x(1)]/(2*area_K);
    Rv_K = eps(K)*(bi*bi'+ci*ci')*area_K;
    Rv(nodes,nodes) = Rv(nodes,nodes)+Rv_K;
end

end