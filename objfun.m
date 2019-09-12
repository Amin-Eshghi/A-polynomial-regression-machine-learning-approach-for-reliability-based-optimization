
function [f,g]=objfun(x) % obj function calculation
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

global cost

W =x(1);
T=x(2);
Sy=21000;
Y=550;
Z=300;
E=15.22*10^6;
f=W*T;
g=[T W];

cost=f;
%fid=fopen('HistConst.txt','a');
%fprintf(fid,['x=[' num2str(x') '], ']);
%fprintf(fid,['f=[' num2str(f) '] ']);
%fprintf(fid,'\n');
%fclose(fid);

return
end