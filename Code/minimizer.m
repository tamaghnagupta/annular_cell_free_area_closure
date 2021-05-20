clc
close all
clear all

A=[];                       % left-hand side of inequality constraint
b=[];                       % right-hand side of inequality constraint
Aeq = [];                   % left-hand side of equality constraint
beq = [];                   % right-hand side of equality constraint
nvars = 4;                  % number of optimization variables

D0_min = 0;                 % lower bound of D0 [mm^2 hr^(-1)]
D0_max = 0.004167;          % upper bound of D0 [mm^2 hr^(-1)]    
D1_min = 0;                 % lower bound of D1 [mm^2 hr^(-1)]
D1_max = 4.166667;          % upper bound of D1 [mm^2 hr^(-1)]
m_min = 0;                  % lower bound of m  [dimensionless]
m_max = 4;                  % upper bound of m  [dimensionless]  
lambda_min = 0.0125;        % lower bound of lambda [hr^(-1)]
lambda_max = 0.0125;        % upper bound of lambda [hr^(-1)]

lb = [D0_min, D1_min, m_min, lambda_min];              % lower bound vector
ub = [D0_max, D1_max, m_max, lambda_max];              % upper bound vector

FitnessFunction = @cost2;   % objective function 
nlcon = @nonlcon;           % nonlinear constraint function

options = gaoptimset('display','iter','TolFun',1e-9,'Generations',500,'PlotFcns',@gaplotpareto);
optimalParameters = gamultiobj(FitnessFunction,nvars,A,b,Aeq,beq,lb,ub,nlcon,options);