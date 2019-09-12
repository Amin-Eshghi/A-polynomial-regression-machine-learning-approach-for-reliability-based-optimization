function [GT_RS, Sen_GT_RS] = ML_RS(u,cid,rt)

% This is a compact code for reliability analysis under uncertainty using a Polynomial Regression Machine Learning approach. 
% The code implements a stochastic response surface method (SRSM) which quantifies the uncertainty in a performance function
% for the purpose of reliability-based design optimization (RBDO). The method first
% uses the Polynomial Regression to approximate the limit state functions (LSF). Then it adds some weights and implements 
% Moving Least Square method to locate the new design point closer to the estimated LSF.
% The two-step process results in a stochastic response surface, with which
% we can apply the Monte Carlo simulation (MCS) to obtain the full probabilistic
% characteristics (e.g., statistical moments, reliability, PDF and quantile) of
% the performance function. One mathematical example is solved
% for demonstration purpose.
%
% VARIABLE DEFINITION
% u: mean vector of random variables; s: standard deviation vector
% cid: constraint number; rt: target reliability
% nv: number of random variables; ns: number of MCS samples
% uniComp: polynomial component function values
% Response_RS: random response values generated using SRSM
% GT_RS: quantile of response at target reliability (CDF) level (used
% as reliability constraint in RBDO)
% Sen_GT_RS: sensivities of quantile w.r.t. means of design variables
%
% FUNCTION DEFINITION
% ML_RS(): main function
% ML_sampling(): identify locations of polynomial regression samples
% FindResponse(): define performance functions

% References:
% 1. Eshghi, Amin Toghi, and Soobum Lee. "Adaptive improved response surface method for reliability-based 
% design optimization." Engineering Optimization (2019): 1-19.



if nargin < 3, rt = []; end
if nargin < 2, cid = []; end
if nargin < 1, u = []; end
if isempty(rt), rt = 0.9987; end        %% Reliability target
if isempty(cid), cid = 2; end           %% Constraint number
if isempty(u), u = [0.57  1.83];end     %% Mean vector of random variables
s = 0.05*[2.5   2.5]';                  %% Standard deviation vector

%===========  Generate MCS Samples based on Uniform Random Seeds =========%
nv = length(u);                         %% Number of input random variables
type = [1 1];                           %% Distribution type: 1: normal; 3:uniform; 6: weibull; 7: lognormal
load UNIFSEED UNIFSEED
ns = size(UNIFSEED,2);
xs = zeros(nv,ns);                      %% Initialization of MCS sample vector
for k = 1:nv
    if type(k)==1                       %% normal distribution
        xs(k,:) = norminv(UNIFSEED(k,:),u(k),s(k));
    end
end

%=====================  Obtain polynomial regression Samples  =======================%
[output] = ML_sampling(u,s,cid);

%================  Response Surface Method (RSM)  ==================%
% Step 1: obtain regression coefficient components
load Xi
x1=Xi(:,1);
x2=Xi(:,2);

Yrsm = output(1,:);
Yrsm=Yrsm';

% Define weights based on distance to the LSF
Weight_Yrsm=Yrsm(:,1)-Yrsm(3,1);
Norm_Weight_Yrsm=Weight_Yrsm./sum(abs(Weight_Yrsm(:,1)));
Size_Norm=max(size(Norm_Weight_Yrsm));
Weight=zeros(Size_Norm,Size_Norm);


for i=1:Size_Norm
   Weight(i,i)=Norm_Weight_Yrsm(i,1);
   if Weight(i,i)==0
      Weight(i,i)=1;
            
   end
end


Xc=[ones(size(x1)),x1,x2,x1.*x1,x2.*x2,x1.*x2];
Xt=Xc';
Z=Xt*Weight*Xc;
V=Xt*Weight*Yrsm;
%regression coefficients
Coeff=Z\V;            

uniComp = zeros(1,ns);
%Obtaining function values using polynomial regression 
for i=1:1:max(size(xs))
   uniComp(1,i)=Coeff(1)+xs(1,i)*Coeff(2)+xs(2,i)*Coeff(3)+Coeff(4)*xs(1,i)*xs(1,i)+Coeff(5)*xs(2,i)*xs(2,i)+Coeff(6)*xs(1,i)*xs(2,i);  
end
Response_RS =uniComp;

%===============  Compute True Responses by Direct MCS  ==================%
%Response_True = FindResponse(xs,cid);

%==================  Conduct Reliability Analysis  =======================%
GT_RS = quantile(Response_RS, rt);    % Quantile of response at target reliability (CDF) level

%==================  Conduct Sensitivity Analysis  =======================%
Sen_GT_RS = zeros(nv,1);
Sen_GUtoU_RS = zeros(nv,1);
Sen_GStoU_RS = zeros(nv,1);
% Compute First-Order Sensitivity of Reliability using Score Function
% Refer to "Stochastic sensitivity analysis by dimensional decomposition and
% score functions", Probabilistic Engineering Mechanics, Volume 24,
% Issue 3, July 2009, Pages 278¨C287
for k = 1:nv
    % Step 1: compute dmu(G)/du = E[G*s], s: score function
    Sen_GUtoU_RS(k) = mean(Response_RS.*((xs(k,:) - u(k))/s(k)^2));
    
    % Step 2: compute dstd(G)/du = (E[G^2*s] - 2*E[G]*E[G*s])/(2*std(G)), s: score function
    Sen_GStoU_RS(k) = mean(Response_RS.^2.*(xs(k,:) - u(k))/s(k)^2) ....
        - 2*mean(Response_RS)*Sen_GUtoU_RS(k);
    Sen_GStoU_RS(k) = Sen_GStoU_RS(k)/2/std(Response_RS);
    
    % Step 3: compute dGT/du =~ dmu(G)/du + 3*dstd(G)/du
    % (an approximation for rt = 0.9978)
    Sen_GT_RS(k) = Sen_GUtoU_RS(k) + 3*Sen_GStoU_RS(k);
end

%==================  Obtain regression Samples  ==========================%

function [output] = ML_sampling(u,s,cid)
global nv m leng Xi
  
nv = length(u);                 %% Dimension of the problem
%m = length(u_loc);             %% Number of samples along each dimension

%Saturated Design (SD)
u_loc_ccd=[-1 0;1 0; 0 0; 0 -1;0 1];
   
leng=max(size(u_loc_ccd));
  
Xi=[];                          % sampling matrix for polynomial regression    



for kk = 1:leng
    
        
    % Identify sample location for RSM
    
    input_ccd(kk,:) = u + (u_loc_ccd(kk,:)'.*s)';
    
    % Get Response values
    xx =  input_ccd(kk,:);
    
            
        Xi=[Xi;xx];
        save Xi Xi
        output(1,kk) = FindResponse(xx,cid);
       
    end
    
%====================== Define Performance Functions ======================%

function [Response]=FindResponse(x,cid)
global c_total count nc nv m leng
count = count + 1;
if isvector(x) == 1
    if cid == 1
        [c ceq] = constraint(x);
        c_total = [c_total;c];
    else          
        c = c_total(count-(cid-1)*leng,:);  %(m+(nv-1)*(m-1)) number of x vectors
    end
    
    if cid == 1
        Response=c(1);
        
        %Response = constraint(x);
    elseif cid == 2
        Response=c(2);
    
    end
    
    if count==2*leng
        count=0;
        c_total=[];
    end
    
end

