%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE MODEL AND THE STAT DISTRIBIUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc
clear all
close all

% PARAMETERS, 1 MONTH IS ONE UNIT

    param0.delta   = 0.0075;
    param0.muH     = 4;
    param0.sigma   = 20;
    param0.b       = 0; % 0;
    param0.gamX    = 0.25;    
    param0.rhoX    = 0.96;
    param0.stdX    = 0.3;% 0.3;
    param0.f       = 0.15;
    param0.psi     = 0.3;
    param0.eta     = 1.3;
    param0.p0      = 0.3;    % mean of p0, coming from unemployment
    param0.p0_sigma= 0.1;    % std of p0, coming from unemployment
    param0.HLboth  = 0;      % HLboth = 1 if switches can happen both from H to L and vice versa; it is 0 if only from H to L
    
    param0.sigma_grid = [17 19 22];

for loop_param  = 1:3
        % CHANGING PARAMETERS
        param = param0;        
        mm_mat = fieldnames( param0 );
        %eval('param.(cell2mat(mm_mat(loop_param)))= param10.(cell2mat(mm_mat(loop_param)))');
        
        % OTHER PARAMETERS
        param.sigma = param0.sigma_grid(loop_param);
        
        param.rho     = 0.004;
        param.muL     = -param.muH;
        param.II      = 5 ;
        param.p1      = param.p0;    % mean of p0, coming from employment
        param.p1_sigma= param.p0_sigma;    % std of p0, coming from employment
        param.Q       = param.gamX*ones(param.II,param.II).*(eye(param.II)==0);
        param.s       = (param.muH-param.muL)/param.sigma;
        [param.Omega param.AA]  = rouwenhorst_approx(param.rhoX,param.stdX,param.II);
        [VV,DD]     = eig(param.Omega');
        xx          = find(diag(DD)<10^(-10));
        param.Omega_bar   = (VV(:,xx)/sum(VV(:,xx)))';  % stationary distribution
        

        
        clear xx VV DD

        % PARAMETERS FOR SIMULATIONS
        param.NN = 10000;
        param.FF = 200;
        param.TT = 15000;
        param.alpha = 0.5;
        param.Nper = 30;
    
    % SOLVE THE SURPLUS
        solve_surplus_bertrand
        seval = ['save sol_betrand_03_' int2str(loop_param) '.mat param sol'];
        eval(seval)

    % STATIONARY DISTRIBUTION
        solve_stat_distribution_p0distr_bertrand
        seval = ['save sol_stat_betrand_03_' int2str(loop_param) '.mat param sol'];
        eval(seval);

    % SIMULATIONS
        simulations_p0distr_bertrand
        seval = ['save simul_betrand_03_' int2str(loop_param) '.mat param sol simul'];
        eval(seval);

    % OUTPUT
        %summarize_output_bertrand
        
end

