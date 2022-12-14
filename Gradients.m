%%% Copied from the book
%%% "The Finite Element Method: Theory, Implementation, and
%%% Practice%: Larson, Bengtzon 
function [area,b,c] = Gradients(x,y)
    area=polyarea(x,y);
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end