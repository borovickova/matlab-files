% SIMULATIONS

% clear all
% close all
% clc

times = zeros(20,1);

%load simul_vars14.mat
% load sol_stat_betrand_DC_1.mat
% load sol_stat_bertrand_DC.mat
% 
%     param.NN = 10000;
%     param.FF = 200;
%     param.TT = 15000;
%     param.alpha = 0.5;
%     param.Nper = 30;
% GET THE PARAMETERS
    mm_mat = fieldnames( param );
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= param.(cell2mat(mm_mat(i)));']);
    end

    mm_mat = fieldnames(sol_DC);
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= sol_DC.(cell2mat(mm_mat(i)));']);
    end

    p1 = p0;
    p1_sigma = p0_sigma;
% distribution of initial beliefs
    fp0 = create_fp0distr(p,p0,p0_sigma);
    fp1 = create_fp0distr(p,p1,p1_sigma);
     
% OTHER VARIABLES
     Nper = 30;
    dt = 1/Nper;
    PP          = expm(Omega*dt);
    csOmega     = [zeros(II,1) cumsum(PP,2)];  
    
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
 
% OLD VERSION    
%     searchers  = NN*(unemrate+psi*(1-unemrate));
%     temp1 = 1/(searchers*f)*Omega_bar./vbar.*(1/eta*Jbar).^(1/(eta-1));
%     qq = (temp1(end)*FF)^(-(eta-1)/eta);    % vacancy-filling rate
%     theta = f/qq;           % vacancy-unemployment ratio
%     VV_theory = searchers*theta;
%     BB = f/theta^(1-alpha); % efficiency in the matching function
%     
%     % number of posted vacancies of firms with productivity A
%     vA = (qq/(cc*eta)*Jbar).^(1/(eta-1));

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
    offer_treshold      = zeros(NN,1);
    endog_sep_treshold  = zeros(NN,1);   
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
    
    mean_belief_p  = zeros(TT,1);
    muH_count = zeros(TT,1);
    muL_count = zeros(TT,1);

    
%     vbar_time = zeros(TT,II);
%     VV_time = zeros(TT,1);
%     f_time = zeros(TT,1);

% INITIAL STAGE
    Ne = floor(NN*(1-unemrate));
    
    firm_name(1:Ne,1)       = randsample(firms,Ne,1);
    firm_name(Ne+1:NN,1)    = 0;
    belief_p(1:Ne,1)        = random_draw_fp0(p,fp0,Ne);
    belief_p(Ne+1:NN,1)     = 0;
    tenure(1:NN,1)          = 0;
    value_mu(1:Ne,1)        = (rand(Ne,1)<belief_p(1:Ne,1))*(muH-muL)+muL;
    value_mu(Ne+1:NN,1)     = 0;
    
    csOmega_bar = [0 cumsum(Omega_bar)];
    %firm_shocks = sum(repmat(rand(1,FF)',1,II+1)>repmat(csOmega_bar,FF,1),2)'; 
    firm_shocks   = randsample([1:1:II], FF, true,vbar) ;
    matrix_shocks = repmat([1:1:II],FF,1);

% save numbers     
    vbar_theory =vbar;
    f_theory = f;

    % EMPIRICIAL DISTRIBUTION
    if p(1)==0
        ptemp = p;
    else
        ptemp = [0 p];
    end

    if p(end)<1
        ptemp = [ptemp 1];
    end

    emp_distr = zeros(Np,II);
    Npp = length(ptemp);
%==========================================================================
% TT1 throw away simulations  
% LOOP

TT1 = 1500;
for tt = 1:TT1
    
    timesind = 1;
    
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
      %fu_endog_hist(tt)=fu_new;
      
      vbar = vv_temp.*vA/VV;
      csvbar = [0 cumsum(vbar)];
      csvbar(end)=1;
      


%0       
% SEPS DUE TO FIRM-LEVEL SHOCKS: those who are under the threshold before
% anything else happens

    % before updating their beliefs
    endog_sep_threshold              = zeros(NN,1);        
    endog_sep_threshold(employed)    = plow(firm_shocks(firm_name(employed)))'; 
    sep_endog_firmshock              = (belief_p < endog_sep_threshold);
    sep_endog_firmshock              = sep_endog_firmshock>0;
    mu_has_changed_sep               = mu_has_changed(sep_endog_firmshock);
    mu_has_changed_sep               = mu_has_changed_sep>0;
    
%1    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % UPDATE TENURE AND BELIEFS
    tenure(employed)   = tenure(employed) + 1;
    dW                 = randn(Ne,1)*sigma*sqrt(dt);
    dZ                 = (dW + value_mu(employed)*dt - belief_p(employed)*muH*dt - (1-belief_p(employed))*muL*dt)/sigma;
    belief_p(employed) = belief_p(employed) + s*belief_p(employed).*(1-belief_p(employed)).*dZ;
    belief_p = min(1,max(0,belief_p));
    
%2    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % SEPARATIONS AND NEW OFFERS
    potential_offer           = sum(repmat(rand(NN,1),1,II+1)>repmat(csvbar,NN,1),2); 
    receive_offer_rand        = rand(NN,1);
    receive_offer             = zeros(NN,1);
    receive_offer(employed)   = receive_offer_rand(employed)  < psi*f*dt;
    receive_offer(unemployed) = receive_offer_rand(unemployed)< f*dt;
    [receive_p0(employed) position_p0(employed)]     = random_draw_fp0(p,fp1,sum(employed));
    [receive_p0(unemployed) position_p0(unemployed)] = random_draw_fp0(p,fp0,sum(unemployed));
    value_mu_new = (rand(Ne+Nu,1)<receive_p0)*(muH-muL)+muL;
    
      
    offer_threshold            = zeros(NN,1);
    accept_offer               = zeros(NN,1)<1;    % this creates a logical array

    current_surplus   = calculate_surplus(p,JJ,belief_p(employed),firm_shocks(firm_name(employed)));
    poaching_surplus  = JJ(position_p0(employed) + Np*(potential_offer(employed)-1));
    poaching_surplus(employed(receive_p0<0)) = 0;   % if p0 is below the lowest possible plow, then poaching surplus is 0
    accept_offer(employed)      = (current_surplus' < poaching_surplus) & receive_offer(employed);
    poaching_surplus_UN  = JJ(position_p0(unemployed) + Np*(potential_offer(unemployed)-1));
    %offer_threshold(unemployed) = plow(potential_offer(unemployed));
    accept_offer(unemployed)    = (poaching_surplus_UN > 0) & receive_offer(unemployed);
    %AO_unem(tt)= sum(accept_offer(unemployed));
        
    sep_exog                    = (rand(NN,1) < delta*dt);
    
    endog_sep_threshold              = zeros(NN,1);        
    endog_sep_threshold(employed)    = plow(firm_shocks(firm_name(employed)))'; 
    sep_endog                        = (belief_p < endog_sep_threshold);
    
%3    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;

    
    % MEASURE SEPARATIONS
    % ORDER: 1) exogenous sep, 2) sep to unemp, 3) sep to to other job

    sep_endog(employed)    = sep_endog(employed).*(1-sep_exog(employed)).*(1-sep_endog_firmshock(employed));    % cannot have endog sep if he had exog sep
    sep_endog(unemployed)  = 0;
    sep_endog = sep_endog  > 0;
    accept_offer(employed) = (1-sep_endog_firmshock(employed)).*(1-sep_exog(employed)).*(1-sep_endog(employed)).*accept_offer(employed);   % cannot accept if endog or exog happened
    accept_offer = accept_offer > 0;
    

%7
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % SEPARATIONS TO UNEMPLOYMENT
    
    %sep_ex_en            = max(max(sep_exog ,sep_endog),sep_endog_firmshock);  
    sep_ex_en            = sep_exog + sep_endog +sep_endog_firmshock;
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
   %change_mu_cum = change_mu_cum + sum(change_mu);
   value_mu(employed(change_mu)) = muL*(value_mu(employed(change_mu))>muL)+ HLboth*muH*(value_mu(employed(change_mu))<muH);
   
   mu_has_changed(employed(change_mu)) = 1;
   
   if tt==100*j
       disp(sprintf('Iter %d emp %d',tt,Ne));
       j = j+1;
   end
   
   firm_shocks = new_shocks;


end

%==========================================================================    
% LOOP

jj = 1;
j  = 1;
max_tenure = 40;

tenure_quar_freq  = zeros(floor(TT/90),max_tenure);
tenure_quar_firms = zeros(floor(TT/90),max_tenure,FF);

%tic; lasttime = 0;
change_mu_cum = 0;

for tt = 1:TT
    
    timesind = 1;
    
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
      %fu_endog_hist(tt)=fu_new;
      
      vbar = vv_temp.*vA/VV;
      csvbar = [0 cumsum(vbar)];
      csvbar(end)=1;
      
%       vbar_time(tt,:) = vbar;
%       VV_time(tt) = VV;
%       f_time(tt) = f;
%       

%0 count workers at the beginning of the period 	
    workers        = (firmsNN == repmat(firm_name,1,FF)); % this identifies which workers belong to each firms
    employment_firm_begin(tt,:) = sum(workers);


%0       
% SEPS DUE TO FIRM-LEVEL SHOCKS: those who are under the threshold before
% anything else happens

    % before updating their beliefs
    endog_sep_threshold              = zeros(NN,1);        
    endog_sep_threshold(employed)    = plow(firm_shocks(firm_name(employed)))'; 
    sep_endog_firmshock              = (belief_p < endog_sep_threshold);
    sep_endog_firmshock              = sep_endog_firmshock>0;
    mu_has_changed_sep               = mu_has_changed(sep_endog_firmshock);
    mu_has_changed_sep               = mu_has_changed_sep>0;
   
    MAT_endog_sep_firmshock = repmat(firms,sum(sep_endog_firmshock),1) == repmat(firm_name(sep_endog_firmshock),1,FF); % this identifies which workers belong to each firms
    endog_sep_firmshocks_firm(tt,:)    = sum(MAT_endog_sep_firmshock,1);
    % WHOSE MU HAS CHANGED?
    if sum(mu_has_changed_sep)>0
        MAT_endog_sep_mu_changed = repmat(firms,sum(mu_has_changed_sep),1) == repmat(firm_name(mu_has_changed_sep),1,FF); % this identifies which workers belong to each firms
        endog_sep_mu_changed_firmshock_firm(tt,:)    = sum(MAT_endog_sep_mu_changed,1);
    end
    
%1    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % UPDATE TENURE AND BELIEFS
    tenure(employed)   = tenure(employed) + 1;
    dW                 = randn(Ne,1)*sigma*sqrt(dt);
    dZ                 = (dW + value_mu(employed)*dt - belief_p(employed)*muH*dt - (1-belief_p(employed))*muL*dt)/sigma;
%     dx_output          = AA(firm_shocks(firm_name(employed)))'*dt + value_mu(employed)*dt + dW;  
    %dZ                = (randn(Ne,1)*sigma*sqrt(dt) + value_mu(employed)*dt - belief_p(employed)*muH*dt - (1-belief_p(employed))*muL*dt)/sigma;
    belief_p(employed) = belief_p(employed) + s*belief_p(employed).*(1-belief_p(employed)).*dZ;
    belief_p = min(1,max(0,belief_p));
    
    % CALCULATE FIRMS OUTPUT
%     MAT_output = repmat(firms,sum(employed),1) == repmat(firm_name(employed),1,FF); % this identifies which workers belong to each firms
%     output_firm(tt,:)    = sum(MAT_output.*repmat(dx_output,1,FF),1);
%     
    mean_belief_p(tt) = mean(belief_p(employed));
    muH_count(tt) = sum(value_mu(employed)==muH);
    muL_count(tt) = sum(value_mu(employed)==muL);
%2    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % SEPARATIONS AND NEW OFFERS
    potential_offer           = sum(repmat(rand(NN,1),1,II+1)>repmat(csvbar,NN,1),2); 
    receive_offer_rand        = rand(NN,1);
    receive_offer             = zeros(NN,1);
    receive_offer(employed)   = receive_offer_rand(employed)  < psi*f*dt;
    receive_offer(unemployed) = receive_offer_rand(unemployed)< f*dt;
    [receive_p0(employed) position_p0(employed)]     = random_draw_fp0(p,fp1,sum(employed));
    [receive_p0(unemployed) position_p0(unemployed)] = random_draw_fp0(p,fp0,sum(unemployed));
    value_mu_new = (rand(Ne+Nu,1)<receive_p0)*(muH-muL)+muL;
    
%     figure(1)
%     subplot(1,2,1)
%     hist(receive_p0(employed),20)
%     subplot(1,2,2)
%     hist(receive_p0(unemployed),20)
%     pause
      
    offer_threshold            = zeros(NN,1);
    accept_offer               = zeros(NN,1)<1;    % this creates a logical array

    current_surplus   = calculate_surplus(p,JJ,belief_p(employed),firm_shocks(firm_name(employed)));
    poaching_surplus  = JJ(position_p0(employed) + Np*(potential_offer(employed)-1));
    poaching_surplus(employed(receive_p0<0)) = 0;   % if p0 is below the lowest possible plow, then poaching surplus is 0
    accept_offer(employed)      = (current_surplus' < poaching_surplus) & receive_offer(employed);
    poaching_surplus_UN  = JJ(position_p0(unemployed) + Np*(potential_offer(unemployed)-1));
    %offer_threshold(unemployed) = plow(potential_offer(unemployed));
    accept_offer(unemployed)    = (poaching_surplus_UN > 0) & receive_offer(unemployed);
    %AO_unem(tt)= sum(accept_offer(unemployed));
        
    sep_exog                    = (rand(NN,1) < delta*dt);
    
    endog_sep_threshold              = zeros(NN,1);        
    endog_sep_threshold(employed)    = plow(firm_shocks(firm_name(employed)))'; 
    sep_endog                        = (belief_p < endog_sep_threshold);
    
%3    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % MEASURE SEPARATIONS
    % ORDER: 1) exogenous sep, 2) sep to unemp, 3) sep to to other job

    sep_endog(employed)    = sep_endog(employed).*(1-sep_exog(employed)).*(1-sep_endog_firmshock(employed));    % cannot have endog sep if he had exog sep
    sep_endog(unemployed)  = 0;
    sep_endog = sep_endog  > 0;
    accept_offer(employed) = (1-sep_endog_firmshock(employed)).*(1-sep_exog(employed)).*(1-sep_endog(employed)).*accept_offer(employed);   % cannot accept if endog or exog happened
    accept_offer = accept_offer > 0;
    
    % ---------------------------------------------------------------------
    
    MAT_endog_sep = repmat(firms,sum(sep_endog),1) == repmat(firm_name(sep_endog),1,FF); % this identifies which workers belong to each firms
    endog_sep_firm(tt,:)    = sum(MAT_endog_sep,1);
    MAT_exog_sep = repmat(firms,sum(sep_exog),1) == repmat(firm_name(sep_exog),1,FF); % this identifies which workers belong to each firms
    exog_sep_firm(tt,:)    = sum(MAT_exog_sep,1);
    MAT_search_sep = repmat(firms,sum(accept_offer),1) == repmat(firm_name(accept_offer),1,FF); % this identifies which workers belong to each firms
    search_sep_firm(tt,:)    = sum(MAT_search_sep,1);
    
    % ---------------------------------------------------------------------
    % WHOSE MU HAS CHANGED
    mu_has_changed(unemployed) = 0;
    sep_endog_mu = (sep_endog.*mu_has_changed)>0;
    sep_exog_mu = (sep_exog.*mu_has_changed)>0;
    accept_offer_mu = (accept_offer.*mu_has_changed)>0;
    
    MAT_endog_sep_mu = repmat(firms,sum(sep_endog_mu),1) == repmat(firm_name(sep_endog_mu),1,FF); % this identifies which workers belong to each firms
    endog_sep_mu_changed_firm(tt,:)    = sum(MAT_endog_sep_mu,1);
    MAT_exog_sep_mu = repmat(firms,sum(sep_exog_mu),1) == repmat(firm_name(sep_exog_mu),1,FF); % this identifies which workers belong to each firms
    exog_sep_mu_changed_firm(tt,:)    = sum(MAT_exog_sep_mu,1);
    MAT_search_sep_mu = repmat(firms,sum(accept_offer_mu),1) == repmat(firm_name(accept_offer_mu),1,FF); % this identifies which workers belong to each firms
    search_sep_mu_changed_firm(tt,:)    = sum(MAT_search_sep_mu,1);
    
%4
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % ---------------------------------------------------------------------

    HIRESall(tt) = sum(accept_offer);
    UEtran(tt)   = sum(accept_offer(unemployed));

%7
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    
    % SEPARATIONS TO UNEMPLOYMENT
    
    %sep_ex_en            = max(max(sep_exog ,sep_endog),sep_endog_firmshock);  
    sep_ex_en            = sep_exog + sep_endog +sep_endog_firmshock;
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
    

      
%9    
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;
    workers        = (firmsNN == repmat(firm_name,1,FF)); % this identifies which workers belong to each firms
    employment_firm_end(tt,:)   = sum(workers);

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
   
   if tt==100*j
       disp(sprintf('Iter %d emp %d changemu %d fu %f8',tt,Ne,change_mu_cum,f));
       j = j+1;
       %EMPIRICAL DISTRIBUTION
          employed          = firm_name > 0;
          Ne = sum(employed);
        Aptemp = repmat(ptemp,Ne,1);
        emp_belief = belief_p(employed);
        belief_temp = repmat(emp_belief,1,Npp);
        p_row = sum(belief_temp>Aptemp,2);
        p_column = firm_shocks(firm_name(employed))';

        new_count = zeros(Np,II);
        new_count(p_row+(p_column-1)*Np) = 1;

        emp_distr = emp_distr+new_count; 
%        
   end
   
   firm_shocks = new_shocks;

%13
%times(timesind) = times(timesind)+toc-lasttime; lasttime = toc; timesind = timesind + 1;

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

simul.mean_belief_p = mean_belief_p;
simul.muH_count = muH_count;
simul.muL_count = muL_count;
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