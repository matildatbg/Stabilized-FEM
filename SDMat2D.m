%%% Copied from the book
%%% "The Finite Element Method: Theory, Implementation, and
%%% Practice%: Larson, Bengtzon 
%%%
%%% Creates streamline diffusion matrix for 2D-problems.

function Sd = SDMat2D(p,t,bx,by)
np=size(p,2);
nt=size(t,2);
Sd=sparse(np,np);
for i=1:nt
loc2glb=t(1:3,i);
x=p(1,loc2glb);
y=p(2,loc2glb);
[area,b,c]=Gradients(x,y);
bxmid=mean(bx(loc2glb));
bymid=mean(by(loc2glb));
SdK=(bxmid*b+bymid*c)*(bxmid*b+bymid*c)'*area;
Sd(loc2glb,loc2glb)=Sd(loc2glb,loc2glb)+SdK;
end