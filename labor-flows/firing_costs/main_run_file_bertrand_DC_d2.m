%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE MODEL AND THE STAT DISTRIBIUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc
clear all
close all

% PARAMETERS, 1 MONTH IS ONE UNIT

% LOAD PARAMETERS FROM THE MODEL WITHOUT DISMISSAL COSTS

%load sol_stat_betrand_bench2.mat
load sol_stat_betrand_09d_2.mat
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
av_value_match  = sum(temp(:))+b/rho;

param.d_grid = av_value_match*[0.005 0.01 0.1];


    % PARAMETERS FOR SIMULATIONS
    param.NN = 10000;
    param.FF = 200;
    param.TT = 15000;
    param.alpha = 0.5;
    param.Nper = 30;
         

for loop_param = 1:3
    param.d = param.d_grid(loop_param);   
    
    % SOLVE THE MODEL
    solve_surplus_bertrand_DC
    seval = ['save sol_betrand_DC_d2' int2str(loop_param) '.mat param sol_DC param_benchmark sol_benchmark'];
    eval(seval)

    % STATIONARY DISTRIBUTION
    solve_stat_distribution_p0distr_bertrand_DC
    seval = ['save sol_stat_betrand_DC_d2' int2str(loop_param) '.mat param sol_DC param_benchmark sol_benchmark'];
    eval(seval);
    
    % SIMULATIONS
    simulations_p0distr_bertrand_DC
    seval = ['save simul_betrand_DC_d2' int2str(loop_param) '.mat param sol_DC param_benchmark sol_benchmark simul'];
    eval(seval);
        
%     for i=1:II
%         temp = sol.JJ(:,i).*(sol.JJ(:,i)<0);
%         [a(i) aa(i)] = max(temp);
%     end

end