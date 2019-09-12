function RBDO_ML()
warning off
% This code consists of a Machine Learning (LM) approach for reliability-based design optimization (RBDO) 
% which uses a polynomial regression technique for reliability analysis under uncertainty. 
% One mathematical example is solved for demonstration purpose. 
% To see more details about this method and other numerical examples you can check out:
%% 
% 1. Eshghi, Amin Toghi, and Soobum Lee. "Adaptive improved response surface method for reliability-based 
% design optimization." Engineering Optimization (2019): 1-19.       
%%                                                                         
% VARIABLE DEFINITION                                                    
% nc: number of constraints ; nv: number of design variables
% rt: target reliability; x0: initial design; xp: previous design
% lb & ub: lower bound and upper bound for design variables               
% iter: variable used to record iteration number of optimization process  
% UNIFSEED: uniform random seeds for reliability analysis using ML
% ns: number of MCS samples
% GT_RS: quantile of response at target reliability (CDF) level (used 
% as reliability constraint in RBDO)
% Sen_GT_RS: sensivities of quantile w.r.t. means of design variables
%                                        
% FUNCTION DEFINITION                                                     
% RBDO_ML(): main function                                                   
% Costfun(): objective function & gradients w.r.t. means of design variables
% frelcon(): Define the probability constraints and gradients             
% ML_RS(): ML method for reliability analysis under uncertainty       
% SHOW(): show the optimization information for every design iteration



global nc rt Iters cost c_total count nv m 
c_total = [];
count = 0;

%===================    Define Optimization Parameters  ===================%
nc = 2; rt = normcdf(3,0,1);
%[Width  Thickness]
x0=[2.5     2.5];

lb=[0  0 ];

ub=[ 5    5];


xp = x0; Iters = 0;
options = optimset('GradConstr','on','GradObj','on','LargeScale','off','Algorithm','sqp');

%==========    Generate Uniform Random Seeds for Reliability Analysis  ========%
nv = length(x0);                     %% Number of input random variables
ns = 1000000;                        %% Number of MCS samples
UNIFSEED = unifrnd(0,1,nv,ns); 
save UNIFSEED UNIFSEED

%=======================    Start Optimization  ==========================%
[x,fval,exitflag,output]=fmincon(@objfun,x0,[],[],[],[],lb,ub,@frelcon,options)

%====================  Define Constraints and Gradiants  =================%
    function [c,ceq,GC,GCeq] = frelcon(x)
        ceq = []; GCeq = []; 
        
        for j = 1:nc
            [GT_RS, Sen_GT_RS] = ML_RS(x,j,rt);
            c(j) = GT_RS;
            GC(:,j) = Sen_GT_RS;
        end
        
        dx = norm(x-xp);
        if  dx > 1d-10  || Iters == 0
            Iters = Iters + 1;
            SHOW(Iters,x,c,GC);
        end
        xp = x;
    end

%===================== Display Iteration Information================%
    function  SHOW(Iters,x,c,GC)
        fprintf(1,'\n********** Iter.%d ***********\n' ,Iters);
        disp(['Des.: ' sprintf('%6.4f  ',x)]);
        disp(['Obj.: ' sprintf('%6.4f',cost)]);
        disp(['Cons.: ' sprintf('%6.4f  ',c)]);
        for k = 1:nv
            if k ==1
                disp(['Sens.: ' sprintf('%6.4f  ',GC(k,:))]);
            else
                disp(['           ' sprintf('%6.4f  ',GC(k,:))]);

            end
        end
        fprintf('\n\n')
    end
%===============================================================%
end
