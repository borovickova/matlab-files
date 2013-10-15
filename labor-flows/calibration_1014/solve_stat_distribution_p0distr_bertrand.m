%function [g ] = solve_stat_distribution_p0distr_bertrand[param,sol]

% calculate stationary distribution for given vbar, f
% clear
% parameters
%load simul_betrand_AD2.mat
% load sol_betrand_12_2.mat

% GET THE PARAMETERS
    mm_mat = fieldnames( param );
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= param.(cell2mat(mm_mat(i)));']);
    end

    mm_mat = fieldnames(sol);
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= sol.(cell2mat(mm_mat(i)));']);
    end
%    
% p0       = 0.3;
% p0_sigma = 0.1;
% p1       = p0;
% p1_sigma = p0_sigma;
% f = 0.2;

JJ = JJ.*(JJ>0);

% % =========================================================================
% % refine grid around 1
% 
% % alternative 1 - geometric;
% leftnode = 980;
% A = p(leftnode)-p(leftnode-1);
% x = 1 - p(leftnode);
% q = 0.96;
% N = log(1-x*(1-q)/(A*q))/log(q);
% refinement = A*q.^[1:N];
% pnew = [p(1:leftnode-1) p(leftnode)+cumsum(refinement) 1];
% 
% % alternative 2 - power
% % leftnode = 980;
% % x0 = p(leftnode);
% % A = p(leftnode)-p(leftnode-1);
% % N = 100;
% % a1 = fminsearch(@(a1) a1search(a1,[x0 N A]), 1);
% % a3 = exp((x0-a1*log(a1))/a1);
% % a2 = a1*a3-a3*x0;
% % [a1 a2 a3]
% % pnew = [p(1:leftnode-1) a1*log(a2+a3*[0:N]*A)]; 
% 
% JJnew = zeros(length(pnew),II);
% for i = 1:II
%     JJnew(:,i) = interp1(p,JJ(:,i),pnew,'spline');
% end
% 
% p = pnew;
% JJ = JJnew;
% Np = length(p);

% =========================================================================

ZZtemp = zeros(Np,II);

s = s*(p<0.9) + max(0,s*(20*(1-p)-1)).*(p>=0.9);

sigmap2 = s.^2.*p.^2.*(1-p).^2;

% differences
dpb = [p(2)-p(1) diff(p)];
dpf = [diff(p) p(end)-p(end-1)];
dpc = [2*p(2)-2*p(1) p(3:end)-p(1:end-2) 2*p(end)-2*p(end-1)];

% distribution of initial beliefs % fp0, fp1 are bin probabilities
fp0 = create_fp0distr(p,p0,p0_sigma);
fp1 = create_fp0distr(p,p1,p1_sigma);

% -------------------------------------------------------------------------
% prepare ODE coefficients
% level term - contributions from 2nd derivative and separation
g_coef = s.^2.*(1-6*p+6*p.^2);
% first-derivative term
gp_coef = s.^2.*(2*p-6*p.^2+4*p.^3);
% second-derivative term
gpp_coef = 1/2*s.^2.*p.^2.*(1-p).^2;

% -------------------------------------------------------------------------
% construct A and b
vbar = ones(1,II)/II;
%vbar = vbar/sum(vbar);

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
                ZZ = (sum((JJ>JJ(k,i)).*repmat(0.5*dpc'.*fp1',1,II))*(vbar')); 
                ZZtemp(k,i) = ZZ - vbar(i)*((0.5*dpc.*fp1)*(JJ(:,i)>JJ(k,i)));

                % above separation threshold
                A(kk,kk) = g_coef(k) - gpp_coef(k)*(1/dpf(k)+1/dpb(k))/dpb(k)...
                           -delta + Omega(i,i)- psi*f*ZZ;
                A(kk,kk-1) = -gp_coef(k)/dpc(k) + gpp_coef(k)/dpb(k)^2;
                A(kk,kk+1) = gp_coef(k)/dpc(k) + gpp_coef(k)/dpf(k)/dpb(k);
                % add switching from different Aj
                for j = 1 : II
                    if ((j ~= i) && (p(k)> HLboth*Q(j,i)) && (p(k)<1-Q(j,i)))
                        porig = (p(k)-HLboth*Q(j,i))/(1-Q(j,i)-HLboth*Q(j,i));
                        porigind = find(p-porig>0,1);
                        if (porigind > 1)
                            wright = (porig - p(porigind-1))/(p(porigind)-p(porigind-1));
                            A(kk,porigind-1+(j-1)*Np) = A(kk,porigind-1+(j-1)*Np) + (1-wright)*Omega(j,i)/(1-Q(j,i)-HLboth*Q(j,i)); %/(1-2*Q(j,i))
                            A(kk,porigind+(j-1)*Np) = A(kk,porigind+(j-1)*Np) + wright*Omega(j,i)/(1-Q(j,i)-HLboth*Q(j,i)); %/(1-2*Q(j,i))
                        end
                    end

                    % add inflow from unemployment
                    A(kk,1+(j-1)*Np:j*Np) = A(kk,1+(j-1)*Np:j*Np) - f*vbar(i)*fp0(k) * dpc*0.5;

                    % add inflow from employment  
                    A(kk,1+(j-1)*Np:j*Np) = A(kk,1+(j-1)*Np:j*Np) + psi*f*vbar(i)*fp1(k) * 0.5*dpc.*(JJ(k,i)>=JJ(:,j)'); % - 1/2*(p==p_bar(j,i)));         
                end

                % RHS
                RHS(kk) = - f*vbar(i)*fp0(k);   %due to inflow from unemployment
            end
        end
    end

    % ---- LAST ROW---------------------------------------------------
%     for i = 1 : II
%         A(i*Np-1,:) = 0;
%         RHS(i*Np-1) = 0;
%         A(i*Np-1,i*Np) = 1;
%     end
    g1is0 = 1;
    if  (g1is0)    
%         for i = 1 : II
%             A(i*Np,:) = 0;
%             RHS(i*Np) = 0;
%             A(i*Np,i*Np) = 1;
%         end
        for i = 1 : II
            A(i*Np,:) = 0;
            A(i*Np,i*Np) = 1;
            A(i*Np,i*Np-1) = -1;
            RHS(i*Np) = -0.0;
        end
    else
        % sum of outflow - sum of inflow = 0
        for i = 1 : II
            rowtofill = i*Np;
            A(rowtofill,:) = 0;
            RHS(rowtofill) = 0;
            % inflows - constant in b, and opposite signs on g
            f0prob = (0.5*dpc.*fp0)*(p>plow(i))';
            RHS(rowtofill) = f*vbar(i)*f0prob;

            for j = 1 : II
                % inflow from unemployment
                A(rowtofill,(j-1)*Np+1:j*Np) =  A(rowtofill,(j-1)*Np+1:j*Np)...
                    + f*vbar(i)*f0prob*dpc*1/2;         

                lower_bound = max(plow(j),(plow(i))/(1-Q(j,i)-HLboth*Q(j,i)));
                %lower_bound = max(plow(i),(plow(j))*(1-Q(j,i)-HLboth*Q(j,i)));
                if (j ~= i)
                    % inflow from employment and switching Aj
%                     offer_accept_in = sum(repmat(fp1'.*dpc'/2,1,Np).*(repmat(JJ(:,j)',Np,1)<repmat(JJ(:,i),1,Np)));
                    %offer_accept_in = sum(repmat(fp1',1,Np).*(repmat(JJ(:,j)',Np,1)<repmat(JJ(:,i),1,Np)));
                    offer_accept_in = (0.5*dpc.*fp1)*(repmat(JJ(:,j)',Np,1)<repmat(JJ(:,i),1,Np));
                    A(rowtofill,(j-1)*Np+1:j*Np) =  A(rowtofill,(j-1)*Np+1:j*Np)...
                        - psi*f*vbar(i)*dpc*1/2.*offer_accept_in...
                        - Omega(j,i)*dpc*1/2.*(p>=lower_bound);
                end
            end
            % outflows
            % outflow to U (exogenous), outflow to E, outflow due to Ai
            A(rowtofill,(i-1)*Np+1:i*Np) =  A(rowtofill,(i-1)*Np+1:i*Np)...
                + delta*dpc*1/2 ...
                + psi*f*dpc*1/2.*ZZtemp(:,i)'...
                - Omega(i,i)*dpc*1/2;
            % outflow to U (endogenous)
            A(rowtofill,(i-1)*Np+plow_pos(i)+1) = A(rowtofill,(i-1)*Np+plow_pos(i)+1) ...
                + 1/2*sigmap2(plow_pos(i))/dpf(plow_pos(i));
            A(rowtofill,(i-1)*Np+plow_pos(i)) = A(rowtofill,(i-1)*Np+plow_pos(i)) ...
                - 1/2*sigmap2(plow_pos(i))/dpf(plow_pos(i));
        end
    end

    % ----------------------------------------------------------
    disp('done')
%     A=sparse(A);
%     RHS = sparse(RHS);
%     g = full(A\RHS);
    g = (A\RHS);
    g = reshape(g',Np,II)';

%     figure
%     plot(p,g)
%     %title('Model with distribution, mean p0, low variance, compicated p=1')
%     title('Model with distribution, mean p0, low variance')

    g1 = g.*(g>0);
    employment  = sum(g*(dpc')/2)
    employment1 = sum(g1*(dpc')/2)
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

            f0prob = (0.5*dpc.*fp0)*(p>plow(i))';

            inf_unem = f*vbar(i)*f0prob*unemrate;
            out_unem = delta*sum(g(i,:).*dpc/2)+ 1/2*sigmap2(plow_pos(i))/dpf(plow_pos(i))*(g(i,plow_pos(i)+1)-g(i,plow_pos(i)));
            out_A = -Omega(i,i)*sum(g(i,:).*dpc/2);

            for j = 1 : II
                % inflow from unemployment
                lower_bound = max(plow(j),(plow(i))/(1-Q(j,i)));

                if (j ~= i)
                    % inflow from employment and switching Aj
                    offer_accept_in = (0.5*dpc.*fp1)*(repmat(JJ(:,j)',Np,1)<repmat(JJ(:,i),1,Np));
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
    hire_from_e = zeros(II,II);
    for i=1:II
        hire_from_u = unemrate/searchers*sum(JJ(:,i)'.*fp0);
        for j=1:II
            T1 = repmat(JJ(:,i),1,Np)-repmat(JJ(:,j)',Np,1);
            T1 = T1.*(T1>0);
            T2 = repmat(fp1',1,Np).*repmat(g(j,:).*dpc/2,Np,1);
            hire_from_e(i,j) = psi*(1-unemrate)/searchers*sum(sum(T1.*T2));
        end

        Jbar(i) = hire_from_u+ sum(hire_from_e(i,:));
    end

    num = Omega_bar.*(Jbar).^(1/(eta-1));
    vbar_new = num/sum(num)
    if max(abs(vbar-vbar_new)) < 10^(-3)
        iter_finish_vbar = 1;
    end
    
    vbar = vbar_new;
end

%----- save variables ------------------------
sol.g = g;
sol.Jbar = Jbar;
sol.vbar = vbar;
sol.unemrate = unemrate;










