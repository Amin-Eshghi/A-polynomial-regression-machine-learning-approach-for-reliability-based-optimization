function [c,ceq]=constraint(x) %Constraint defined here

% Numerical example to demonstrate the performance of the proposed method. This example is provided with detailed explaination at:
% References:
% 1. Eshghi, Amin Toghi, and Soobum Lee. "Adaptive improved response surface method for reliability-based 
% design optimization." Engineering Optimization (2019): 1-19.

% Cantilever beam 
% This example demonstrates a cantilever beam design problem. The design goal
% is to minimize its weight subject to multiple constraints on displacement and stress. 
% The cantilever beam is in vertical and lateral bending which is loaded at the tip by the vertical and lateral loads Y and Z, respectively. 
% Beam length is 100 in. The width w and thickness t of the cross section are deterministic design variables. 
% With the assumption of constant density, the objective function is equivalent to minimizing the cross-sectional area defined as f = w * t. 
% Two nonlinear constraints are considered in the formulation. The first constraint (G1) is related to yielding at the fixed end of the beam; 
% the other constraint (G2) deals with the tip displacement which is allowed to be less than D0=2.2535 in.

L=100;
D0=2.2535;
        
W =x(1);
T=x(2);
Sy=21000;
Y=550;
Z=300;
E=15.22*10^6;

c(1)=-(Sy-((600*Y/(W*T^2))+(600*Z/(W^2*T))));
c(2)=-(D0-(4*L^3*(((Y/(T^2))^2+(Z/(W^2))^2)^0.5)/(E*W*T)));

ceq=[];

fid=fopen('HistConst.txt','a');
fprintf(fid,['x=[' num2str(x) '], ']);
fprintf(fid,['c=[' num2str(c) '] ']);
fprintf(fid,'\n');
fclose(fid);

return
end