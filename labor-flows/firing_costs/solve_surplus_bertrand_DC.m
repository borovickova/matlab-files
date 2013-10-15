%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL MODEL WITH SEARCH ON THE JOB AND SHOCKS TO MATCH QUALITY
% ADDITIVE SHOCKS
% Bertrand wage setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GET THE PARAMETERS
mm_mat = fieldnames( param );
for i = 1 : length(mm_mat)
    eval([cell2mat(mm_mat(i)) '= param.(cell2mat(mm_mat(i)));']);
end

% SOLVE THE SURPLUS
tic;

iter_finish = zeros(1,II);
iter = zeros(1,II);
           
plow_pos = ones(1,II);

disp(sprintf('======================================================='))
    
% OTHER VARIABLES
Np_initial = 1000;

% ITERATE
iter_plow_all_finish = 0;
iter_plow_all = 1;

pmin   = zeros(1,II);
pmax   = ones(1,II); 
plow = pmin;
plow_all_old = ones(1,II);
    
    while(~iter_plow_all_finish)
                
            for iplow_iter = 1:II
                                
                %disp(sprintf('Iterations for shock %d',iplow_iter))    
                
                iter_finish(iplow_iter) = 0;
                iter(iplow_iter) = 1;
                                
%                 if (~iter_finish(iplow_iter))
%                     if (iplow_iter == II)
%                         pmin(iplow_iter) = 0;
%                     else
%                         pmin(iplow_iter) = pmin(iplow_iter+1);
%                     end
%                     if (iplow_iter == 1)
%                         pmax(iplow_iter) = 1;
%                     else
%                         pmax(iplow_iter) = pmax(iplow_iter-1);
%                     end
%                     plow(iplow_iter) = (pmin(iplow_iter) + pmax(iplow_iter))/2;
%                 end

   
                    pmin(iplow_iter) = 0;
                    pmax(iplow_iter) = 1;
%                     if (iplow_iter == 1)
%                         pmax(iplow_iter) = 1;
%                     else
%                         pmax(iplow_iter) = pmax(iplow_iter-1);
%                     end
                    plow(iplow_iter) = (pmin(iplow_iter) + pmax(iplow_iter))/2;
        

                while(~iter_finish(iplow_iter))
                    %disp(sprintf('Plow iteration %d',iter(iplow_iter)))

                    % CREATE GRID
                    p = create_grid_bertrand(plow,Np_initial); 
                    Np = length(p);
                    [plow_pos, JJ, times] =  solve_JJ_loop_bertrand_DC(param,p,Np,plow);
                    
                    %total_times = zeros(10,1);

                    plow_pos(iplow_iter) = min(Np-1,plow_pos(iplow_iter));
                    Dfirst = (JJ(plow_pos(iplow_iter)+1,iplow_iter)-JJ(plow_pos(iplow_iter),iplow_iter))/(p(plow_pos(iplow_iter)+1)-p(plow_pos(iplow_iter)));
                    if (Dfirst > 0)
                        pmax(iplow_iter) = p(plow_pos(iplow_iter));
                    else
                        pmin(iplow_iter) = p(plow_pos(iplow_iter));
                    end 
                   
                    plow(iplow_iter) = (pmin(iplow_iter) + pmax(iplow_iter))/2;

                    if (abs(pmin(iplow_iter)-pmax(iplow_iter)) < 10^(-5)) ||  iter(iplow_iter)>30
                        iter_finish(iplow_iter)= 1;
                    end

                    iter(iplow_iter)  = iter(iplow_iter)+1;
                end
                
                disp(sprintf('Iteration: plow all %d, shock %d, shock plow %d, time %f',iter_plow_all,iplow_iter,iter(iplow_iter), toc));
                disp(plow)
            end 

        % VERIFY SMOOTH-PASTING FOR ALL i

        for i=1:II
            plow_pos(i) = min(Np-5,plow_pos(i));
            Dfirst_all(i) = (JJ(plow_pos(i)+1,i)-JJ(plow_pos(i),i))/(p(plow_pos(i)+1)-p(plow_pos(i)));
        end

        Dfirst_all
        Dfirst_all = Dfirst_all.*(plow>0.0001);
        Dfirst_all = Dfirst_all.*(plow<0.9999);
        Dplow = abs(plow_all_old -plow);
        
%         pause
        if  max(abs(Dfirst_all))<10^(-2)  || max(Dplow)<10^(-5) | (iter_plow_all >15)
            iter_plow_all_finish = 1;   
        end

%         if  max(abs(pmax-pmin))<10^(-6)  | (iter_plow_all >50)
%             iter_plow_all_finish = 1;   
%         end
%         figure(1)
%         plot(p,JJ)
%        pause

        iter_plow_all = iter_plow_all + 1;
        plow_all_old = plow;
    end
    
%----- save variables ------------------------
sol_DC.p = p;
sol_DC.JJ = JJ;
sol_DC.plow = plow;
sol_DC.plow_pos =plow_pos;
sol_DC.Np = Np;

    
