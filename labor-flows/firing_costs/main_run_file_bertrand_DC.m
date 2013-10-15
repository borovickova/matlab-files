%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE MODEL AND THE STAT DISTRIBIUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc
clear all
close all

% PARAMETERS, 1 MONTH IS ONE UNIT

% LOAD PARAMETERS FROM THE MODEL WITHOUT DISMISSAL COSTS

load sol_stat_betrand_12_2.mat
param_benchmark = param;
sol_benchmark = sol;

mm_mat = fieldnames( param );
for i = 1 : length(mm_mat)
    eval([cell2mat(mm_mat(i)) '= param.(cell2mat(mm_mat(i)));']);
end

% mm_mat = fieldnames( sol );
% for i = 1 : length(mm_mat)
%     eval([cell2mat(mm_mat(i)) '= sol.(cell2mat(mm_mat(i)));']);
% end

p = sol_benchmark.p;
dpc = [2*p(2)-2*p(1) p(3:end)-p(1:end-2) 2*p(end)-2*p(end-1)];
temp = (sol_benchmark.JJ').*(repmat(dpc,II,1)/2.*sol_benchmark.g);
av_value_firm  = sum(temp(:));

param.d_grid = av_value_firm*[0.005 0.01 0.1];

%     param0.delta   = 0;
%     param0.muH     = 5;
%     param0.sigma   = 55;
%     param0.b       = 0;
%     param0.gamX    = 0.25;    
%     param0.rhoX    = 0.8;
%     param0.stdX    = 0.3;
%     param0.f       = 0.15;
%     param0.psi     = 0.17;
%     param0.eta     = 2;
%     param0.p0      = 0.6;    % mean of p0, coming from unemployment
%     param0.p0_sigma= 0.1;    % std of p0, coming from unemployment
%     param0.HLboth  = 1;      % HLboth = 1 if switches can happen both from H to L and vice versa; it is 0 if only from H to L
%     
% 
%     param = param0;        
%     % OTHER PARAMETERS
%     param.rho     = 0.004;
%     param.muL     = -param.muH;
%     param.II      = 5 ;
%     param.p1      = param.p0;    % mean of p0, coming from employment
%     param.p1_sigma= param.p0_sigma;    % std of p0, coming from employment
%     param.Q       = param.gamX*ones(param.II,param.II).*(eye(param.II)==0);
%     param.s       = (param.muH-param.muL)/param.sigma;
%     [param.Omega param.AA]  = rouwenhorst_approx(param.rhoX,param.stdX,param.II);
%     [VV,DD]     = eig(param.Omega');
%     xx          = find(diag(DD)<10^(-10));
%     param.Omega_bar   = (VV(:,xx)/sum(VV(:,xx)))';  % stationary distribution
%     clear xx VV DD

    % PARAMETERS FOR SIMULATIONS
    param.NN = 10000;
    param.FF = 200;
    param.TT = 15000;
    param.alpha = 0.5;
    param.Nper = 30;
         

for loop_param = 2:3
%    param.d = param.d_grid(loop_param);   
    
%     % SOLVE THE MODEL
%     solve_surplus_bertrand_DC
%     seval = ['save sol_betrand_DC_x' int2str(loop_param) '.mat param sol_DC param_benchmark sol_benchmark'];
%     eval(seval)
% 
%     % STATIONARY DISTRIBUTION
%     solve_stat_distribution_p0distr_bertrand_DC
%     seval = ['save sol_stat_betrand_DC_x' int2str(loop_param) '.mat param sol_DC param_benchmark sol_benchmark'];
%     eval(seval);
    
    % SIMULATIONS
    seval = ['load sol_stat_betrand_DC_x' int2str(loop_param) '.mat'];
    eval(seval);
    simulations_p0distr_bertrand_DC
    seval = ['save simul_betrand_DC_xx' int2str(loop_param) '.mat param sol_DC param_benchmark sol_benchmark simul'];
    eval(seval);
        
%     for i=1:II
%         temp = sol.JJ(:,i).*(sol.JJ(:,i)<0);
%         [a(i) aa(i)] = max(temp);
%     end

end