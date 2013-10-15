function [plow_pos, JJ, times] =  solve_JJ_loop_bertrand_DC(param,p,Np,plow);

    % retrieve parameters
    mm_mat = fieldnames( param );
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= param.(cell2mat(mm_mat(i)));']);
    end
   
times = zeros(10,1);
% tic;

sigmap2 = s^2*p.^2.*(1-p).^2;

% differences
diffx = p(2:end)-p(1:end-1);
dpb = [p(2)-p(1) diffx];
dpf = [diffx p(end)-p(end-1)];
%dpc = [p(2)-p(1) p(3:end)-p(1:end-2) p(end)-p(end-1)];

% find indices of plow thresholds and p0

% prepare grid matrices
% bigB = sparse(Np*II,Np*II);
% RHS = sparse(Np*II,1);
bigB = eye(Np*II,Np*II);
RHS = zeros(Np*II,1);

% disp(sprintf('1: %d %d',size(bigB)));

% times(1) = toc; tic;

% -------------------------------------------------------------------------
% prepare ODE coefficients
% level term - contributions from 2nd derivative and separation
plow_pos = zeros(1,II);
for i = 1 : II
       temp = p - plow(i);
       temp = temp.*(temp<0);
       [a aa] = max(temp);
       plow_pos(i) = aa;
end


% times(2) = toc; tic;

% -------------------------------------------------------------------------
% construct A and b

% prepare grid matrices


    mup = muH*p + muL*(1-p);
    
for k = 1 : Np-1
    for i = 1 : II
        % matrix row
        kk = k+(i-1)*Np;
        if (k <= plow_pos(i))
            % below separation threshold
            bigB(kk,kk) = 1;
            RHS(kk) = -d;
%             times(3) = times(3) + toc; tic;
        else
            % above separation threshold
             bigB(kk,kk)   = rho+delta -Omega(i,i) + 0.5*sigmap2(k)*(1/dpf(k)+1/dpb(k))/dpf(k);
             bigB(kk,kk+1) = -0.5*sigmap2(k)/(dpf(k)^2);
             bigB(kk,kk-1) = -0.5*sigmap2(k)/dpf(k)/dpb(k);
             RHS(kk) = AA(i)+ mup(k)-b;
%             times(4) = times(4) + toc; tic;
            % add switching from different Aj
            for j = 1 : II
                if (j ~= i)
                     pprime = p(k)*(1-Q(i,j))+ HLboth*(1-p(k))*Q(i,j) ;
                     pprimeind = find(p-pprime>0,1);
%                    times(5) = times(5) + toc; tic;
%                    times(6) = times(6) + toc; tic;
                    if (pprimeind > 1)
                          wright = (pprime - p(pprimeind-1))/(p(pprimeind)-p(pprimeind-1));
%                          times(7) = times(7) + toc; tic;
                         bigB(kk,pprimeind-1+(j-1)*Np) = -(1-wright)*Omega(i,j);
                         bigB(kk,pprimeind+(j-1)*Np) = -wright*Omega(i,j);
                    end
%                     times(8) = times(8) + toc; tic;
                end
%                 times(9) = times(9) + toc; tic;
            end
        end 
    end
end

% disp(sprintf('2: %d %d',size(bigB)));
% last row
k = Np;

    for i = 1 : II
        % matrix row
        kk = k+(i-1)*Np;
            % above separation threshold
             bigB(kk,kk) = rho + delta -Omega(i,i) ;
             RHS(kk) = AA(i)+ mup(k)-b;
            % add switching from different Aj
            for j = 1 : II
                 if (j ~= i)
                     pprime = p(k)*(1-Q(i,j))+ HLboth*(1-p(k))*Q(i,j);
                     pprimeind = find(p-pprime>0,1);
                    if (pprimeind > 1)
                          wright = (pprime - p(pprimeind-1))/(p(pprimeind)-p(pprimeind-1));
                         bigB(kk,pprimeind-1+(j-1)*Np) = -(1-wright)*Omega(i,j);
                         bigB(kk,pprimeind+(j-1)*Np) = -wright*Omega(i,j);
                    end
                end
            end
    end
 
% disp(sprintf('3: %d %d',size(bigB)));
    bigBsparse = sparse(bigB);
    RHSsparse = sparse(RHS);
    JJtemp2 = (bigBsparse\(RHSsparse))';
%     JJtemp2 = (bigB\(RHS))';
    JJtemp = full(JJtemp2);
    JJ     = reshape(JJtemp,Np,II);

% times(10) = toc; tic;

% whos
end