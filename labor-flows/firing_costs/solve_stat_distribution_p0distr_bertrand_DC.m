%function [g ] = solve_stat_distribution_p0distr_bertrand[param,sol]

% calculate stationary distribution for given vbar, f
% clear
% parameters
% load sol_betrand.mat
% load sol_betrand_DC_1.mat

% GET THE PARAMETERS
    mm_mat = fieldnames( param );
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= param.(cell2mat(mm_mat(i)));']);
    end

    mm_mat = fieldnames(sol_DC);
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= sol_DC.(cell2mat(mm_mat(i)));']);
    end
   
% p0      = 0.3;
% p0_sigma= 0.1;
% f = 0.15;

JJ = JJ.*(JJ>0);
ZZtemp = zeros(Np,II);

% distribution of initial beliefs
fp0 = create_fp0distr(p,p0,p0_sigma);
fp1 = create_fp0distr(p,p1,p1_sigma);

sigmap2 = s^2*p.^2.*(1-p).^2;

% differences
dpb = [p(2)-p(1) diff(p)];
dpf = [diff(p) p(end)-p(end-1)];
dpc = [2*p(2)-2*p(1) p(3:end)-p(1:end-2) 2*p(end)-2*p(end-1)];


% -------------------------------------------------------------------------
% prepare ODE coefficients
% level term - contributions from 2nd derivative and separation
g_coef = s^2*(1-6*p+6*p.^2);
% first-derivative term
gp_coef = s^2*(2*p-6*p.^2+4*p.^3);
% second-derivative term
gpp_coef = 1/2*s^2*p.^2.*(1-p).^2;

% -------------------------------------------------------------------------
% construct A and b
vbar = ones(1,II)/II;
%vbar = vbar/sum(vbar);

% -------------------------------------------------------------------------
% RECOVER THE PARAMETER B in the matching function
% Omega_bar_bench, unemrate_bench, Jbar_bench

Omega_bar_bench = param_benchmark.Omega_bar;
unemrate_bench = sol_benchmark.unemrate;
Jbar_bench = sol_benchmark.Jbar;
vbar_bench = sol_benchmark.vbar;
f_bench = param_benchmark.f;
% NEW c1, c0
c1 = eta-1;
searchers  = NN*(unemrate_bench + psi*(1-unemrate_bench));
    temp1 = 1/(searchers*f_bench)*Omega_bar_bench./vbar_bench.*Jbar_bench.^(1/c1);;
    qq = (temp1(end)*FF)^(-c1/(1+c1));    % vacancy-filling rate
    theta = f/qq;           % vacancy-unemployment ratio
    VV_theory = searchers*theta;
    BB = f/theta^(1-alpha); % efficiency in the matching function
    
   
% -------------------------------------------------------------------------


iter_finish_f = 0;

while (~iter_finish_f)
    iter_finish_vbar = 0;
while(~iter_finish_vbar)
    % prepare grid matrices
    A   = zeros(Np*II,Np*II);
    RHS = zeros(Np*II,1);

    for k = 1 : Np-1
        for i = 1 : II
            % matrix row
            kk = k+(i-1)*Np;
            if (k <= plow_pos(i))
                % below separation threshold
                A(kk,kk) = 1;
                RHS(kk) = 0;
            else
                % outflow to other employment
                ZZ = sum(sum((JJ>JJ(k,i)).*repmat(fp1'.*dpc'/2,1,II)).*vbar);
                ZZtemp(k,i) = ZZ - vbar(i)*sum((JJ(:,i)>JJ(k,i)).*fp1'.*dpc'/2);

                % above separation threshold
                A(kk,kk) = g_coef(k) - gpp_coef(k)*(1/dpf(k)+1/dpb(k))/dpb(k)...
                           -delta + Omega(i,i)- psi*f*ZZ;
                A(kk,kk-1) = -gp_coef(k)/dpc(k) + gpp_coef(k)/dpb(k)^2;
                A(kk,kk+1) = gp_coef(k)/dpc(k) + gpp_coef(k)/dpf(k)/dpb(k);
                % add switching from different Aj
                for j = 1 : II
                    if ((j ~= i) && (p(k)> HLboth*Q(j,i)) && (p(k)<1-Q(j,i)))
                        pprime = (p(k)-HLboth*Q(j,i))/(1-Q(j,i)-HLboth*Q(j,i));
                        pprimeind = find(p-pprime>0,1);
                        if (pprimeind > 1)
                            wright = (pprime - p(pprimeind-1))/(p(pprimeind)-p(pprimeind-1));
                            A(kk,pprimeind-1+(j-1)*Np) = (1-wright)*Omega(j,i)/(1-Q(j,i)-HLboth*Q(j,i)); %/(1-2*Q(j,i))
                            A(kk,pprimeind+(j-1)*Np) = wright*Omega(j,i)/(1-Q(j,i)-HLboth*Q(j,i)); %/(1-2*Q(j,i))
                        end
                    end

                    % add inflow from unemployment
                    A(kk,1+(j-1)*Np:j*Np) = A(kk,1+(j-1)*Np:j*Np) - f*vbar(i)*fp0(k)*dpc*0.5.*(JJ(:,j)'>=0);

                    % add inflow from employment  
                    A(kk,1+(j-1)*Np:j*Np) = A(kk,1+(j-1)*Np:j*Np) + psi*f*vbar(i)*fp1(k)*dpc*1/2.*(JJ(k,i)>=JJ(:,j)'); % - 1/2*(p==p_bar(j,i)));         
                end

                % RHS
                RHS(kk) = RHS(kk) - f*vbar(i)*fp0(k);   %due to inflow from unemployment
            end
        end
    end

    % ---- LAST ROW---------------------------------------------------
    g1is0 = 1;
    if  g1is0==1    
        for i = 1 : II
            A(i*Np,:) = 0;
            RHS(i*Np) = 0;
            A(i*Np,i*Np)=1;
        end
    else
        % sum of inflow = sum of outflow
        for i = 1 : II
            % inflows - constant in b, and opposite signs on g
            f0prob = sum(fp0.*(p>plow(i)).*dpc/2);
            RHS(i*Np) = f*vbar(i)*f0prob;

            for j = 1 : II
                % inflow from unemployment
                lower_bound = max(plow(j),(plow(i))/(1-Q(j,i)-HLboth*Q(j,i)));
                %lower_bound = max(plow(j),(plow(i))/(1-Q(j,i)-HLboth*Q(j,i)));
                A(i*Np,(j-1)*Np+1:j*Np) =  A(i*Np,(j-1)*Np+1:j*Np)...
                    + f*vbar(i)*f0prob*dpc*1/2;         
                if (j ~= i)
                    % inflow from employment and switching Aj
                    offer_accept_in = sum(repmat(fp1'.*dpc'/2,1,Np).*(repmat(JJ(:,j)',Np,1)<repmat(JJ(:,i),1,Np)));
                    A(i*Np,(j-1)*Np+1:j*Np) =  A(i*Np,(j-1)*Np+1:j*Np)...
                        - psi*f*vbar(i)*dpc*1/2.*offer_accept_in...
                        - Omega(j,i)*dpc*1/2.*(p>lower_bound);
                end
            end
            % outflows
            % outflow to U (exogenous), outflow to E, outflow due to Ai
            A(i*Np,(i-1)*Np+1:i*Np) =  A(i*Np,(i-1)*Np+1:i*Np)...
                + delta*dpc*1/2 ...
                + psi*f*dpc*1/2.*ZZtemp(:,i)'...
                - Omega(i,i)*dpc*1/2;
            % outflow to U (endogenous)
            A(i*Np,(i-1)*Np+plow_pos(i)+1) = A(i*Np,(i-1)*Np+plow_pos(i)+1) ...
                + 1/2*sigmap2(plow_pos(i))/dpf(plow_pos(i));
            A(i*Np,(i-1)*Np+plow_pos(i)) = A(i*Np,(i-1)*Np+plow_pos(i)) ...
                - 1/2*sigmap2(plow_pos(i))/dpf(plow_pos(i));
        end
    end

    % ----------------------------------------------------------
    disp('done')
    A=sparse(A);
    RHS = sparse(RHS);
    g = full(A\RHS);
    g = reshape(g',Np,II)';

%     figure
%     plot(p,g)
%     %title('Model with distribution, mean p0, low variance, compicated p=1')
%     title('Model with distribution, mean p0, low variance')

    g1 = g.*(g>0);
    employment = sum(sum(repmat(dpc/2,II,1).*g))
    employment1 = sum(sum(repmat(dpc/2,II,1).*g1))
    % g1 = g;
    % g1(:,1007:1010) = 0;
    % employment1 = sum(sum(repmat(dpc/2,II,1).*g1))

    %---- check whether inflow(i)==outlow(i)-------------------------------

    unemrate = 1 - employment
    %unemrate = 0.15;

        for i = 1 : II
            % inflows - constant in b, and opposite signs on g
            inf_emp = zeros(1,II);
            inf_A = zeros(1,II);
            out_emp= zeros(1,II);

            f0prob = sum(fp0.*(p>plow(i)).*dpc/2);

            inf_unem = f*vbar(i)*f0prob*unemrate;
            out_unem = delta*sum(g(i,:).*dpc/2)+ 1/2*sigmap2(plow_pos(i))/dpf(plow_pos(i))*(g(i,plow_pos(i)+1)-g(i,plow_pos(i)));
            out_A = -Omega(i,i)*sum(g(i,:).*dpc/2);

            for j = 1 : II
                % inflow from unemployment
                lower_bound = max(plow(j),(plow(i))/(1-Q(j,i)-HLboth*Q(j,i)));

                if (j ~= i)
                    % inflow from employment and switching Aj
                    offer_accept_in = sum(repmat(fp1'.*dpc'/2,1,Np).*(repmat(JJ(:,j)',Np,1)<repmat(JJ(:,i),1,Np)));
                    inf_emp(j) = psi*f*vbar(i)*sum(g(j,:).*offer_accept_in.*dpc/2);
                    inf_A(j) = sum(Omega(j,i)*g(j,:).*dpc/2.*(p>lower_bound));
                    out_emp(j) = psi*f*vbar(j)*sum(g(i,:).*dpc/2.*ZZtemp(:,i)');

                end
            end
            flows(i) = inf_unem+ sum(inf_emp)+ sum(inf_A)-sum(out_emp)-out_A-out_unem;
        end
        flows


    %-----verify guess on vbar ------------------------------------------------
    % I need to calculate firm's expected value of a new hire

    searchers = unemrate+ psi*(1-unemrate);
    Jbar = zeros(1,II);
    for i=1:II
        hire_from_u = unemrate/searchers*sum(JJ(:,i)'.*fp0.*dpc/2);
        T1 = repmat(JJ(:,i),1,Np)-repmat(JJ(:,j)',Np,1);
        T1 = T1.*(T1>0);
        T2 = repmat((fp1.*dpc/2)',1,Np).*repmat(g(j,:).*dpc/2,Np,1);
        hire_from_e = psi*(1-unemrate)/searchers*sum(sum(T1.*T2));

        Jbar(i) = hire_from_u+hire_from_e;
    end

    num = Omega_bar.*(Jbar).^(1/(eta-1));
    vbar_new = num/sum(num)
    if max(abs(vbar-vbar_new)) < 10^(-3)
        iter_finish_vbar = 1;
    end
    
    vbar = vbar_new;
end

    % CHECK IF GUESS FOR f WAS CORRECT
    
    searchers  = NN*(unemrate+psi*(1-unemrate));
    temp1 =  FF/searchers*Omega_bar./vbar.*Jbar.^(1/c1)*BB^(-1/alpha);
    qq = temp1(end)^(-1/(1/alpha + 1/c1));    % vacancy-filling rate
    
    f_new = qq^(1-1/alpha)*BB^(1/alpha);
    
    
    if abs(f-f_new) < 10^(-3)
       iter_finish_f = 1;
    end
    
    f = f_new;
end

vA = (qq*Jbar).^(1/c1);

%----- save variables ------------------------
sol_DC.g = g;
sol_DC.Jbar = Jbar;
sol_DC.vbar = vbar;
sol_DC.unemrate = unemrate;
sol_DC.vA = vA;
sol_DC.f = f_new;










