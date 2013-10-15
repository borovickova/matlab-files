% SIMULATIONS

% clear all
% close all
% clc
% 
% tic;
% 
% times = zeros(20,1);

% load sol_stat_bertrand_12_2.mat 

% GET THE PARAMETERS
    mm_mat = fieldnames( param );
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= param.(cell2mat(mm_mat(i)));']);
    end

    mm_mat = fieldnames(sol);
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= sol.(cell2mat(mm_mat(i)));']);
    end
    
% p0       = 0.3;
% p0_sigma = 0.1;
% p1       = p0;
% p1_sigma = p0_sigma;

% distribution of initial beliefs
%     fp0 = create_fp0distr(p,p0,p0_sigma);
%     fp1 = create_fp0distr(p,p1,p1_sigma);
     
% OTHER VARIABLES
%     Nper = 30;
    dt = 1/Nper;
    PP          = expm(Omega*dt);
    csOmega     = [zeros(II,1) cumsum(PP,2)];  
    Tdrop = 0;
    s_const =  (muH-muL)/sigma;
% FIRMS AND PEOPLE
%     NN = 10000;
%     FF = 200;
%     NN = 50;
%     FF = 10;
    firms = [1:1:FF];
    firmsNN = repmat(firms,NN,1);
    
% RECOVER PARAMETERS NEEDED FOR SIMULATIONS
% ?? IS THIS PART CORRECT FOR THE BERTRAND AND P0 ??
    % matching function alpha; BB
    % vacancy posting function
    % alpha = 0.5;    % elasticity of the matching function
    c0 = 1;         % normalize vacancy posting constant to 1
    c1 = eta-1;

%    OLD: c,eta
%     searchers  = NN*(unemrate+psi*(1-unemrate));
%     temp1 = 1/(searchers*f)*Omega_bar./vbar.*(1/eta*Jbar).^(1/(eta-1));
%     qq = (temp1(end)*FF)^(-(eta-1)/eta);    % vacancy-filling rate
%     theta = f/qq;           % vacancy-unemployment ratio
%     VV_theory = searchers*theta;
%     BB = f/theta^(1-alpha); % efficiency in the matching function
    
    % number of posted vacancies of firms with productivity A
%    vA = (qq/(eta)*Jbar).^(1/(eta-1));
 
% NEW c1, c0
    searchers  = NN*(unemrate+psi*(1-unemrate));
    temp1 = 1/(searchers*f)*Omega_bar./vbar.*Jbar.^(1/c1);
    qq = (temp1(end)*FF)^(-c1/(1+c1));    % vacancy-filling rate
    theta = f/qq;           % vacancy-unemployment ratio
    VV_theory = searchers*theta;
    BB = f/theta^(1-alpha); % efficiency in the matching function

    vA = (Jbar*qq/c0).^(1/c1);

    % PERIODS
    % TT  = 100000;
    % TT  = 50000;
    % TT = 5000;
    
% DEFINE VARIABLES
    offer_threshold     = zeros(NN,1);
    endog_sep_threshold = zeros(NN,1);   
    receive_p0          = zeros(NN,1);
    position_p0         = zeros(NN,1);
    mu_has_changed      = zeros(NN,1);

% TIME VARIABLES    
    employed_time       = zeros(TT,1);
    exog_sep_firm       = zeros(TT,FF);
    endog_sep_firm      = zeros(TT,FF);
    search_sep_firm     = zeros(TT,FF);
    workers_hired_firm  = zeros(TT,FF);
    employment_firm_begin = zeros(TT,FF);
    employment_firm_end = zeros(TT,FF);
    output_firm   = zeros(TT,FF);
    HIRESall              = zeros(TT,1);
    UEtran              = zeros(TT,1);
    fu_endog_hist       = zeros(TT,1);  
    endog_sep_mu_changed_firm = zeros(TT,FF);
    exog_sep_mu_changed_firm = zeros(TT,FF);
    search_sep_mu_changed_firm = zeros(TT,FF);
    endog_sep_firmshocks_firm = zeros(TT,FF);
    endog_sep_mu_changed_firmshock_firm = zeros(TT,FF);
     emp_distr = zeros(Np,II);
%     vbar_time = zeros(TT,II);
%     VV_time = zeros(TT,1);
%     f_time = zeros(TT,1);

% INITIAL STAGE
    Ne = floor(NN*(1-unemrate));
    firm_name               = zeros(NN,1);
    firm_name(1:Ne,1)       = randsample(firms,Ne,1);
    belief_p                = zeros(NN,1);
%     belief_p(1:Ne,1)        = random_draw_fp0(p,fp0,Ne);
    belief_p(1:Ne,1)        = random_draw_fp0_normal(Ne,p0,p0_sigma);
    tenure(1:NN,1)          = 0;
    value_mu                = zeros(NN,1);
    value_mu(1:Ne,1)        = (rand(Ne,1)<belief_p(1:Ne,1))*(muH-muL)+muL;
    
    csOmega_bar = [0 cumsum(Omega_bar)];
    %firm_shocks = sum(repmat(rand(1,FF)',1,II+1)>repmat(csOmega_bar,FF,1),2)'; 
    firm_shocks   = randsample([1:1:II], FF, true,Omega_bar) ;
    matrix_shocks = repmat([1:1:II],FF,1);

% save numbers     
    vbar_theory =vbar;

    f_theory = f;
% pause

%==========================================================================    
% variables for empirical distribution
if p(1)==0
    ptemp = p;
else
    ptemp = [0 p];
end

if p(end)<1
    ptemp = [ptemp 1];
end

Npp = length(ptemp);

%==========================================================================  
% OBSERVATIONS TO DROP
TT1 = 2000;
for tt = 1:TT1
    
    timesind = 1;

    firm_name_old = firm_name;
      
    % IDENTIFY EMPLOYED AND UNEMPLOYED
    employed    = (firm_name >  0);
    unemployed  = (firm_name == 0);
    Ne  = sum(employed);
    Nu  = sum(unemployed);
    employed_time(tt) = Ne;
    
    % CALCULATE ENDOGENOUS VARIABLES
    vv_temp =  sum(matrix_shocks==repmat(firm_shocks',1,II));
    VV = sum(vv_temp.*vA);
    searchers = Nu + psi*Ne;
    f = BB*(VV/searchers)^(1-alpha);
    % fu_endog_hist(tt)=fu_new;

    % distribution of vacancies
    vbar = vv_temp.*vA/VV;
    csvbar = [0 cumsum(vbar)];
    csvbar(end)=1;

%1    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % UPDATE TENURE AND BELIEFS
    tenure(employed)   = tenure(employed) + 1;
    dW                 = randn(Ne,1)*sigma*sqrt(dt);
    emp_belief_p       = belief_p(employed);
    dZ                 = (dW + value_mu(employed)*dt - emp_belief_p*muH*dt - (1-emp_belief_p)*muL*dt)/sigma;
    dx_output          = AA(firm_shocks(firm_name(employed)))'*dt + value_mu(employed)*dt + dW;  
    %dZ                = (randn(Ne,1)*sigma*sqrt(dt) + value_mu(employed)*dt - belief_p(employed)*muH*dt - (1-belief_p(employed))*muL*dt)/sigma;
    % ZRIHLYTTTTTTTTTTTTTTTTTTTTTTT zmikunde
    s                   = s_const*(emp_belief_p<0.9) + max(0,s_const*(20*(1-emp_belief_p)-1)).*(emp_belief_p>=0.9);
    belief_p(employed)  = emp_belief_p + s.*emp_belief_p.*(1-emp_belief_p).*dZ;
    belief_p            = min(1,max(0,belief_p));
    
    
%2    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % SEPARATIONS AND NEW OFFERS
    potential_offer           = sum(repmat(rand(NN,1),1,II+1)>repmat(csvbar,NN,1),2); 
    receive_offer_rand        = rand(NN,1);
    receive_offer             = zeros(NN,1);
    receive_offer(employed)   = receive_offer_rand(employed)  < psi*f*dt;
    receive_offer(unemployed) = receive_offer_rand(unemployed)< f*dt;
    [receive_p0(employed)]    = random_draw_fp0_normal(sum(employed),p1,p1_sigma);
    [receive_p0(unemployed)]  = random_draw_fp0_normal(sum(unemployed),p0,p0_sigma);
    value_mu_new = (rand(Ne+Nu,1)<receive_p0)*(muH-muL)+muL;
    
      
    offer_threshold            = zeros(NN,1);
    accept_offer               = zeros(NN,1)<1;    % this creates a logical array

    current_surplus     = calculate_surplus(p,JJ,belief_p(employed),firm_shocks(firm_name(employed)));
%     poaching_surplus  = JJ(position_p0(employed) + Np*(potential_offer(employed)-1));
    poaching_surplus    = calculate_surplus(p,JJ,receive_p0(employed),potential_offer(employed));
%     poaching_surplus(employed(receive_p0<0)) = 0;   % if p0 is below the lowest possible plow, then poaching surplus is 0
    accept_offer(employed)      = (current_surplus < poaching_surplus) & (receive_offer(employed)');
    offer_threshold(unemployed) = plow(potential_offer(unemployed));
    accept_offer(unemployed)    = (receive_p0(unemployed) >= offer_threshold(unemployed)) & receive_offer(unemployed);
    %AO_unem(tt)= sum(accept_offer(unemployed));
        
    sep_exog                         = (rand(NN,1) < delta*dt);
    endog_sep_threshold              = zeros(NN,1);        
    endog_sep_threshold(employed)    = plow(firm_shocks(firm_name(employed)))'; 
    sep_endog                        = (belief_p < endog_sep_threshold);

%3    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % MEASURE SEPARATIONS
    % ORDER: 1) exogenous sep, 2) sep to unemp, 3) sep to to other job

    sep_endog(employed)    = sep_endog(employed).*(1-sep_exog(employed));    % cannot have endog sep if he had exog sep
    sep_endog(unemployed)  = 0;
    sep_endog = sep_endog  > 0;
    accept_offer(employed) = (1-sep_exog(employed)).*(1-sep_endog(employed)).*accept_offer(employed);   % cannot accept if endog or exog happened
    accept_offer = accept_offer > 0;
    
%7
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % SEPARATIONS TO UNEMPLOYMENT
    
    %sep_ex_en            = max(max(sep_exog ,sep_endog),sep_endog_firmshock);  
    sep_ex_en            = sep_exog + sep_endog;
    sep_ex_en            = sep_ex_en>0;
    belief_p(sep_ex_en)  = 0;
    tenure(sep_ex_en)    = 0;
    firm_name(sep_ex_en) = 0;
    value_mu(sep_ex_en)  = 0; 
    mu_has_changed_sep(sep_ex_en) = 0;
    
    % HIRING: people who find job assign to firms
    

    move_to_shock = accept_offer.*potential_offer;
%8
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    for ishock = 1:II
        Nworker = sum(move_to_shock ==ishock);
        Nfirmspool = firms(firm_shocks ==ishock);
%       pause
        
        move_to_firms = randsample(Nfirmspool,Nworker,1);
        firm_name(move_to_shock == ishock) = move_to_firms;
        workers_hired = repmat(firms,Nworker,1)== repmat(firm_name(move_to_shock == ishock),1,FF);
        workers_hired_firm(tt,:) = workers_hired_firm(tt,:) +  sum(workers_hired,1);
        
        % set new variables
        tenure(move_to_shock == ishock)     = 0;
        belief_p(move_to_shock ==ishock)    = receive_p0(move_to_shock ==ishock);
        value_mu(move_to_shock ==ishock)    = value_mu_new(move_to_shock ==ishock);
        mu_has_changed_sep(move_to_shock ==ishock) = 0;
    end    

%11
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
   
   % NEW FIRM SHOCKS

   xx              = rand(1,FF);
   new_shocks      = sum(repmat(xx',1,II+1)>csOmega(firm_shocks,:),2)';          
   change_shock    = new_shocks~=firm_shocks;
   
   employed          = firm_name > 0;
   Ne                = sum(employed);
   gammas            = (Q((new_shocks(firm_name(employed))-1)*size(Q,1) + firm_shocks(firm_name(employed))) )';   % these are individual gammas for given A(t)-A(t+1)
   belief_p(employed)= (1-gammas).* belief_p(employed) + HLboth*gammas.*(1-belief_p(employed)); % update belief; if A does not change, gamma = 0

%12
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
   
   change_mu = (rand(Ne,1)<gammas);
   change_temp = HLboth*ones(Ne,1) + (1-HLboth)*(value_mu(employed)>muL); % if HLboth=0, then consider only changes if mu==muH
   change_mu = change_mu & (change_temp>0);
   %othermu   = muL*(value_mu(employed)>muL)+ muH*(value_mu(employed)<muH);
   value_mu(employed(change_mu)) = muL*(value_mu(employed(change_mu))>muL)+ HLboth*muH*(value_mu(employed(change_mu))<muH);
   
   mu_has_changed(employed(change_mu)) = 1;
   
   firm_shocks = new_shocks;
   
% 13 
% firm-shock separations    

   endog_sep_threshold              = zeros(NN,1);        
   endog_sep_threshold(employed)    = plow(firm_shocks(firm_name(employed)))'; 
   sep_endog_firmshock              = (belief_p < endog_sep_threshold);
   sep_endog_firmshock              = sep_endog_firmshock>0;
   mu_has_changed_sep               = mu_has_changed(sep_endog_firmshock);
   mu_has_changed_sep               = mu_has_changed_sep>0;
  
   % SEPARATE THE WORKERS
   belief_p(sep_endog_firmshock)  = 0;
   tenure(sep_endog_firmshock)    = 0;
   firm_name(sep_endog_firmshock) = 0;
   value_mu(sep_endog_firmshock)  = 0; 
   mu_has_changed_sep(sep_endog_firmshock) = 0;
    
   employed = firm_name>0;
   unemployed = firm_name==0;
   Ne = sum(employed);

    if (mod(tt,100) == 0)
       disp(sprintf('Iter %d emp %d fu %f8; %f secs.',tt,Ne,f, toc));       
    end
  

end

%==========================================================================  

% LOOP
jj = 1;
max_tenure = 40;

tenure_quar_freq  = zeros(floor(TT/90),max_tenure);
tenure_quar_firms = zeros(floor(TT/90),max_tenure,FF);

%tic; lasttime = 0;
change_mu_cum = 0;

for tt = 1:TT
    
    timesind = 1;

    firm_name_old = firm_name;
      
    % IDENTIFY EMPLOYED AND UNEMPLOYED
    employed    = (firm_name >  0);
    unemployed  = (firm_name == 0);
    Ne  = sum(employed);
    Nu  = sum(unemployed);
    employed_time(tt) = Ne;
    
    % CALCULATE ENDOGENOUS VARIABLES
    vv_temp =  sum(matrix_shocks==repmat(firm_shocks',1,II));
    VV = sum(vv_temp.*vA);
    searchers = Nu + psi*Ne;
    f = BB*(VV/searchers)^(1-alpha);
    % fu_endog_hist(tt)=fu_new;

    % distribution of vacancies
    vbar = vv_temp.*vA/VV;
    csvbar = [0 cumsum(vbar)];
    csvbar(end)=1;
      
%0 count workers at the beginning of the period 	
%     workers        = (firmsNN == repmat(firm_name,1,FF)); % this identifies which workers belong to each firms
    workerstab = tabulate(firm_name);
%     employment_firm_begin(tt,:) = sum(workers);
    employment_firm_begin(tt,1:size(workerstab,1)) = workerstab(:,2)';

%1    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % UPDATE TENURE AND BELIEFS
    tenure(employed)   = tenure(employed) + 1;
    dW                 = randn(Ne,1)*sigma*sqrt(dt);
    emp_belief_p       = belief_p(employed);
    dZ                 = (dW + value_mu(employed)*dt - emp_belief_p*muH*dt - (1-emp_belief_p)*muL*dt)/sigma;
    dx_output          = AA(firm_shocks(firm_name(employed)))'*dt + value_mu(employed)*dt + dW;  
    %dZ                = (randn(Ne,1)*sigma*sqrt(dt) + value_mu(employed)*dt - belief_p(employed)*muH*dt - (1-belief_p(employed))*muL*dt)/sigma;
    % ZRIHLYTTTTTTTTTTTTTTTTTTTTTTT zmikunde
    s                   = s_const*(emp_belief_p<0.9) + max(0,s_const*(20*(1-emp_belief_p)-1)).*(emp_belief_p>=0.9);
    belief_p(employed)  = emp_belief_p + s.*emp_belief_p.*(1-emp_belief_p).*dZ;
    belief_p            = min(1,max(0,belief_p));
    
    % CALCULATE FIRMS OUTPUT
    MAT_output = repmat(firms,sum(employed),1) == repmat(firm_name(employed),1,FF); % this identifies which workers belong to each firms
    output_firm(tt,:)    = sum(MAT_output.*repmat(dx_output,1,FF),1);

%     MAT_output = tabulate(firm_name(sep_endog_firmshock));
%     if (size(MAT_endog_sep_firmshock,1) > 0)
%         endog_sep_firmshocks_firm(tt,1:size(MAT_endog_sep_firmshock,1)) = MAT_endog_sep_firmshock(:,2)';
%     end

    
%2    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % SEPARATIONS AND NEW OFFERS
    potential_offer           = sum(repmat(rand(NN,1),1,II+1)>repmat(csvbar,NN,1),2); 
    receive_offer_rand        = rand(NN,1);
    receive_offer             = zeros(NN,1);
    receive_offer(employed)   = receive_offer_rand(employed)  < psi*f*dt;
    receive_offer(unemployed) = receive_offer_rand(unemployed)< f*dt;
    [receive_p0(employed)]    = random_draw_fp0_normal(sum(employed),p1,p1_sigma);
    [receive_p0(unemployed)]  = random_draw_fp0_normal(sum(unemployed),p0,p0_sigma);
    value_mu_new = (rand(Ne+Nu,1)<receive_p0)*(muH-muL)+muL;
    
%     figure(1)
%     subplot(1,2,1)
%     hist(receive_p0(employed),20)
%     subplot(1,2,2)
%     hist(receive_p0(unemployed),20)
%     pause
      
    offer_threshold            = zeros(NN,1);
    accept_offer               = zeros(NN,1)<1;    % this creates a logical array

    current_surplus     = calculate_surplus(p,JJ,belief_p(employed),firm_shocks(firm_name(employed)));
%     poaching_surplus  = JJ(position_p0(employed) + Np*(potential_offer(employed)-1));
    poaching_surplus    = calculate_surplus(p,JJ,receive_p0(employed),potential_offer(employed));
%     poaching_surplus(employed(receive_p0<0)) = 0;   % if p0 is below the lowest possible plow, then poaching surplus is 0
    accept_offer(employed)      = (current_surplus < poaching_surplus) & (receive_offer(employed)');
    offer_threshold(unemployed) = plow(potential_offer(unemployed));
    accept_offer(unemployed)    = (receive_p0(unemployed) >= offer_threshold(unemployed)) & receive_offer(unemployed);
    %AO_unem(tt)= sum(accept_offer(unemployed));
        
    sep_exog                         = (rand(NN,1) < delta*dt);
    endog_sep_threshold              = zeros(NN,1);        
    endog_sep_threshold(employed)    = plow(firm_shocks(firm_name(employed)))'; 
    sep_endog                        = (belief_p < endog_sep_threshold);

%3    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % MEASURE SEPARATIONS
    % ORDER: 1) exogenous sep, 2) sep to unemp, 3) sep to to other job

    sep_endog(employed)    = sep_endog(employed).*(1-sep_exog(employed));    % cannot have endog sep if he had exog sep
    sep_endog(unemployed)  = 0;
    sep_endog = sep_endog  > 0;
    accept_offer(employed) = (1-sep_exog(employed)).*(1-sep_endog(employed)).*accept_offer(employed);   % cannot accept if endog or exog happened
    accept_offer = accept_offer > 0;
    
    % ---------------------------------------------------------------------
    
%     MAT_endog_sep = repmat(firms,sum(sep_endog),1) == repmat(firm_name(sep_endog),1,FF); % this identifies which workers belong to each firms
%     endog_sep_firm(tt,:)    = sum(MAT_endog_sep,1);
%     MAT_exog_sep = repmat(firms,sum(sep_exog),1) == repmat(firm_name(sep_exog),1,FF); % this identifies which workers belong to each firms
%     exog_sep_firm(tt,:)    = sum(MAT_exog_sep,1);
%     MAT_search_sep = repmat(firms,sum(accept_offer),1) == repmat(firm_name(accept_offer),1,FF); % this identifies which workers belong to each firms
%     search_sep_firm(tt,:)    = sum(MAT_search_sep,1);
% ---------------------------------------------------------------------
%     % WHOSE MU HAS CHANGED
%     mu_has_changed(unemployed)  = 0;
%     sep_endog_mu                = (sep_endog.*mu_has_changed)>0;
%     sep_exog_mu                 = (sep_exog.*mu_has_changed)>0;
%     accept_offer_mu             = (accept_offer.*mu_has_changed)>0;
    
%     MAT_endog_sep_mu = repmat(firms,sum(sep_endog_mu),1) == repmat(firm_name(sep_endog_mu),1,FF); % this identifies which workers belong to each firms
%     endog_sep_mu_changed_firm(tt,:)    = sum(MAT_endog_sep_mu,1);
%     MAT_exog_sep_mu = repmat(firms,sum(sep_exog_mu),1) == repmat(firm_name(sep_exog_mu),1,FF); % this identifies which workers belong to each firms
%     exog_sep_mu_changed_firm(tt,:)    = sum(MAT_exog_sep_mu,1);
%     MAT_search_sep_mu = repmat(firms,sum(accept_offer_mu),1) == repmat(firm_name(accept_offer_mu),1,FF); % this identifies which workers belong to each firms
%     search_sep_mu_changed_firm(tt,:)    = sum(MAT_search_sep_mu,1);
    
    % WHOSE MU HAS CHANGED
    mu_has_changed(unemployed)  = 0;
    sep_endog_mu                = (sep_endog.*mu_has_changed)>0;
    sep_exog_mu                 = (sep_exog.*mu_has_changed)>0;
    accept_offer_mu             = (accept_offer.*mu_has_changed)>0;
    
    if sum(sep_endog)>0
        MAT_endog_sep = tabulate(firm_name(sep_endog));
        endog_sep_firm(tt,1:size(MAT_endog_sep,1)) = MAT_endog_sep(:,2);
        if sum(sep_endog_mu)>0
            MAT_endog_sep_mu = tabulate(firm_name(sep_endog_mu));
            endog_sep_mu_changed_firm(tt,1:size(MAT_endog_sep_mu,1)) = MAT_endog_sep_mu(:,2);
        end
    end
    if sum(sep_exog)>0
        MAT_exog_sep = tabulate(firm_name(sep_exog));
        exog_sep_firm(tt,1:size(MAT_exog_sep,1))   = MAT_exog_sep(:,2);
        if sum(sep_exog_mu)>0
            MAT_exog_sep_mu = tabulate(firm_name(sep_exog_mu));
            exog_sep_mu_changed_firm(tt,1:size(MAT_exog_sep_mu,1)) = MAT_exog_sep_mu(:,2);
        end
    end
    if sum(accept_offer(employed))>0
        MAT_search_sep = tabulate(firm_name(accept_offer(employed)));
        search_sep_firm(tt,1:size(MAT_search_sep,1)) = MAT_search_sep(:,2);
        if sum(accept_offer_mu(employed))>0
            MAT_search_sep_mu = tabulate(firm_name(accept_offer_mu(employed)));
            search_sep_mu_changed_firm(tt,1:size(MAT_search_sep_mu,1)) = MAT_search_sep_mu(:,2);
        end
    end
    % ---------------------------------------------------------------------
     
%4
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % ---------------------------------------------------------------------

    HIRESall(tt) = sum(accept_offer);
    UEtran(tt)   = sum(accept_offer(unemployed));

%7
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % SEPARATIONS TO UNEMPLOYMENT
    
    %sep_ex_en            = max(max(sep_exog ,sep_endog),sep_endog_firmshock);  
    sep_ex_en            = sep_exog + sep_endog;
    sep_ex_en            = sep_ex_en>0;
    belief_p(sep_ex_en)  = 0;
    tenure(sep_ex_en)    = 0;
    firm_name(sep_ex_en) = 0;
    value_mu(sep_ex_en)  = 0; 
    mu_has_changed_sep(sep_ex_en) = 0;
    
    % HIRING: people who find job assign to firms
    
    workers_hired_firm(tt,:) = zeros(1,FF);  
    move_to_shock = accept_offer.*potential_offer;
%8
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    for ishock = 1:II
        Nworker = sum(move_to_shock ==ishock);
        Nfirmspool = firms(firm_shocks ==ishock);
%       pause
        
        move_to_firms = randsample(Nfirmspool,Nworker,1);
        firm_name(move_to_shock == ishock) = move_to_firms;
        workers_hired = repmat(firms,Nworker,1)== repmat(firm_name(move_to_shock == ishock),1,FF);
        workers_hired_firm(tt,:) = workers_hired_firm(tt,:) +  sum(workers_hired,1);
        
        % set new variables
        tenure(move_to_shock == ishock)     = 0;
        belief_p(move_to_shock ==ishock)    = receive_p0(move_to_shock ==ishock);
        value_mu(move_to_shock ==ishock)    = value_mu_new(move_to_shock ==ishock);
        mu_has_changed_sep(move_to_shock ==ishock) = 0;
    end
    

%========================================================================    
%     disp(sprintf('----------------'))
% 
%     disp(sprintf('SEP_EN_FS SEP_EN ACC_OFFER SEP_EX FIRM_NAME'))
%     [sep_endog_firmshock sep_endog accept_offer sep_exog firm_name]
%     disp(sprintf('EMP begin, %d %d', employment_firm_begin(tt,:)))
%     disp(sprintf('HIRES,     %d %d',workers_hired_firm(tt,:)))
%     disp(sprintf('SEP EN FS, %d %d',endog_sep_firmshocks_firm(tt,:)))
%     disp(sprintf('SEP EN,    %d %d',endog_sep_firm(tt,:)))
%     disp(sprintf('SEP EX,    %d %d',exog_sep_firm(tt,:)))
%     disp(sprintf('SEARCH,    %d %d',search_sep_firm(tt,:)))
%     disp(sprintf('EMP end,   %d %d',employment_firm_end(tt,:)))
%     pause
%     
%========================================================================   
    
%10
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    if (mod(tt,90) == 0)
        temp4 = floor((tenure+89.5)/90);% ((tenure>90*(kk-1)) & (tenure<=90*kk));
        temp5 = repmat(temp4,1,FF).*workers;
        for kk = 1:max_tenure
            tenure_quar_freq(jj,kk) = sum(temp4==kk);  
            tenure_quar_firms(jj,kk,:)= sum(temp5==kk);
        end
        jj = jj+1;
    end

%11
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
   
   % NEW FIRM SHOCKS

   xx              = rand(1,FF);
   new_shocks      = sum(repmat(xx',1,II+1)>csOmega(firm_shocks,:),2)';          
   change_shock    = new_shocks~=firm_shocks;
   
   employed          = firm_name > 0;
   Ne                = sum(employed);
   gammas            = (Q((new_shocks(firm_name(employed))-1)*size(Q,1) + firm_shocks(firm_name(employed))) )';   % these are individual gammas for given A(t)-A(t+1)
   belief_p(employed)= (1-gammas).* belief_p(employed) + HLboth*gammas.*(1-belief_p(employed)); % update belief; if A does not change, gamma = 0

%12
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
   
   change_mu = (rand(Ne,1)<gammas);
   change_temp = HLboth*ones(Ne,1) + (1-HLboth)*(value_mu(employed)>muL); % if HLboth=0, then consider only changes if mu==muH
   change_mu = change_mu & (change_temp>0);
   %othermu   = muL*(value_mu(employed)>muL)+ muH*(value_mu(employed)<muH);
   change_mu_cum = change_mu_cum + sum(change_mu);
   value_mu(employed(change_mu)) = muL*(value_mu(employed(change_mu))>muL)+ HLboth*muH*(value_mu(employed(change_mu))<muH);
   
   mu_has_changed(employed(change_mu)) = 1;
   
   firm_shocks = new_shocks;
   
% 13 
% firm-shock separations    

   endog_sep_threshold              = zeros(NN,1);        
   endog_sep_threshold(employed)    = plow(firm_shocks(firm_name(employed)))'; 
   sep_endog_firmshock              = (belief_p < endog_sep_threshold);
   sep_endog_firmshock              = sep_endog_firmshock>0;
   mu_has_changed_sep               = mu_has_changed(sep_endog_firmshock);
   mu_has_changed_sep               = mu_has_changed_sep>0;
   
   MAT_endog_sep_firmshock = tabulate(firm_name(sep_endog_firmshock));
   if (size(MAT_endog_sep_firmshock,1) > 0)
        endog_sep_firmshocks_firm(tt,1:size(MAT_endog_sep_firmshock,1)) = MAT_endog_sep_firmshock(:,2)';
   end
    % WHOSE MU HAS CHANGED?
   if sum(mu_has_changed_sep)>0
        MAT_endog_sep_mu_changed = tabulate(firm_name(mu_has_changed_sep));
        endog_sep_mu_changed_firmshock_firm(tt,1:size(MAT_endog_sep_mu_changed,1)) = MAT_endog_sep_mu_changed(:,2)';
   end
  
   % SEPARATE THE WORKERS
   belief_p(sep_endog_firmshock)  = 0;
   tenure(sep_endog_firmshock)    = 0;
   firm_name(sep_endog_firmshock) = 0;
   value_mu(sep_endog_firmshock)  = 0; 
   mu_has_changed_sep(sep_endog_firmshock) = 0;
    
   employed = firm_name>0;
   unemployed = firm_name==0;
   Ne = sum(employed);

    if (mod(tt,100) == 0)

        xx_fired = ((firm_name_old > 0) & (firm_name == 0));
        xx_hired = ((firm_name_old == 0) & (firm_name > 0));
        xx_switched = ((abs(firm_name_old - firm_name) > 0) & (~xx_hired) & (~xx_fired));
        xx_firmtypes = histc(firm_shocks(firm_name(employed)),[-0.5:1:II+0.5]);
       
        xx_empbelief = belief_p(employed);
        xx_belief_firm1 = xx_empbelief(firm_shocks(firm_name(employed)) == 1);
        xx_belief_firm1 = sort(xx_belief_firm1,'ascend');
        xx_belief_firm2 = xx_empbelief(firm_shocks(firm_name(employed)) == 2);
        xx_belief_firm2 = sort(xx_belief_firm2,'ascend');
        xx_belief_firm3 = xx_empbelief(firm_shocks(firm_name(employed)) == 3);
        xx_belief_firm3 = sort(xx_belief_firm3,'ascend');
        xx_belief_firm4 = xx_empbelief(firm_shocks(firm_name(employed)) == 4);
        xx_belief_firm4 = sort(xx_belief_firm4,'ascend');
        xx_belief_firm5 = xx_empbelief(firm_shocks(firm_name(employed)) == 5);
        xx_belief_firm5 = sort(xx_belief_firm5,'ascend');
%    belief_p
   
       
        disp(sprintf('Iter %d emp %d changemu %d fu %f8 fired %d hired %d switch %d; %f secs.',tt,Ne,change_mu_cum,f, sum(xx_fired), sum(xx_hired),sum(xx_switched), toc));
       
%        hist(firm_shocks)
%        pause
    end
    
%9    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    workers        = (firmsNN == repmat(firm_name,1,FF)); % this identifies which workers belong to each firms
    employment_firm_end(tt,:)   = sum(workers);

    
    % empirical distribution of beliefs
    employed = firm_name>0;
    Ne = sum(employed);

    %if (tt*dt > Tdrop)
        emp_belief = belief_p(employed);
        for i = 1 : II
            curdist = histc(emp_belief(firm_shocks(firm_name(employed)) == i),ptemp);
            emp_distr(:,i) = emp_distr(:,i)+curdist(2:end); 
        end
    %end

       
    %        figure(fighist);
    %        xx_firm = histc(firm_name,[-0.5:1:FF+0.5]);
    %        subplot(3,3,1);
    %        bar([0.5:1:FF-0.5],xx_firm(2:end-1),'histc');
    %        subplot(3,3,2);
    %        plot(sort(xx_firm(2:end-1),'ascend'));
    %        subplot(3,3,3);
    %        bar([0.5:1:II-0.5],xx_firmtypes(2:end-1),'histc');
    %       subplot(3,3,4); hold off;
    %        plot(xx_belief_firm1); hold on; plot([1 length(xx_belief_firm1)],[plow(1) plow(1)],'r:');
    %        subplot(3,3,5); hold off;
    %        plot(xx_belief_firm2); hold on; plot([1 length(xx_belief_firm2)],[plow(2) plow(2)],'r:');
    %        subplot(3,3,6); hold off;
    %        plot(xx_belief_firm3); hold on; plot([1 length(xx_belief_firm3)],[plow(3) plow(3)],'r:');
    %        subplot(3,3,7); hold off;
    %        plot(xx_belief_firm4); hold on; plot([1 length(xx_belief_firm4)],[plow(4) plow(4)],'r:');
    %        subplot(3,3,8); hold off;
    %        plot(xx_belief_firm5); hold on; plot([1 length(xx_belief_firm5)],[plow(5) plow(5)],'r:');
    %        pause


    %13
    %times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    %firm_shocks = new_shocks;
end

% ik = 5;
% figure
% hold on
% plot(vbar_time(:,ik))
% plot(ones(1,TT)*sol.vbar(ik),'r-')
% plot(ones(1,TT)*mean(vbar_time(:,ik)),'g-')
    

%==========================================================================
% SAVE VARIABLES

% simul.Nper = Nper;
 simul.dt = dt;
% simul.NN = NN;
% simul.FF = FF;
% simul.TT = TT;
simul.vbar = vbar;

simul.tenure_quar_freq  = tenure_quar_freq;
simul.tenure_quar_firms = tenure_quar_firms;

simul.employed_time       = employed_time;
simul.exog_sep_firm       = exog_sep_firm;
simul.endog_sep_firm      = endog_sep_firm;
simul.search_sep_firm     = search_sep_firm;
simul.workers_hired_firm  = workers_hired_firm;
simul.employment_firm_begin = employment_firm_begin;
simul.employment_firm_end = employment_firm_end;
simul.output_firm   = output_firm;
simul.HIRESall            = HIRESall;
simul.UEtran              = UEtran;
simul.fu_endog_hist       = fu_endog_hist;
simul.workers_hired_firm  = workers_hired_firm;

simul.endog_sep_mu_changed_firm = endog_sep_mu_changed_firm;
simul.exog_sep_mu_changed_firm = exog_sep_mu_changed_firm;
simul.search_sep_mu_changed_firm = search_sep_mu_changed_firm;
simul.endog_sep_firmshocks_firm = endog_sep_firmshocks_firm;
simul.endog_sep_mu_changed_firmshock_firm = endog_sep_mu_changed_firmshock_firm;

simul.emp_distr = emp_distr;

%==========================================================================
% CHECK IF EMPLOYMENT IS CONSISTENT WITH HIRES AND SEPARATIONS
% x1 =  employment_firm_end(:) - employment_firm_begin(:);
% x2 = workers_hired_firm(:) - (exog_sep_firm(:)+ endog_sep_firm(:)+search_sep_firm(:)+endog_sep_firmshocks_firm(:));
% sum(x1~=x2)
% 
% x3 = HIRESall;
% x4 = sum(workers_hired_firm,2);
% sum(x3~=x4)

%dpc = [2*p(2)-2*p(1) p(3:end)-p(1:end-2) 2*p(end)-2*p(end-1)];

%dpc_x = repmat(dpc'/2,1,II);
%empirical_dist = emp_distr/sum(emp_distr(:))./dpc_x;
%figure
%plot(p,empirical_dist)
%hold on
%plot(p,g')