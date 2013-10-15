% SOLVE WORKER'S SURPLUS, GIVEN THE SURPLUS OF THE MATCH

clear all
close all
clc

load sol_stat_bertrand.mat
 
% GET THE PARAMETERS
    mm_mat = fieldnames( param );
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= param.(cell2mat(mm_mat(i)));']);
    end

    mm_mat = fieldnames(sol);
    for i = 1 : length(mm_mat)
        eval([cell2mat(mm_mat(i)) '= sol.(cell2mat(mm_mat(i)));']);
    end
    
p1 = p0;
p1_sigma = p0_sigma;

fp0 = create_fp0distr(p,p0,p0_sigma);
fp1 = create_fp0distr(p,p1,p1_sigma);

% OTHER USEFUL VALUES    

sigmap2 = s^2*p.^2.*(1-p).^2;

% differences
diffx = p(2:end)-p(1:end-1);
dpb = [p(2)-p(1) diffx];
dpf = [diffx p(end)-p(end-1)];
dpc = [p(2)-p(1) p(3:end)-p(1:end-2) p(end)-p(end-1)]; 

% FIX A VALUE w
% GUESS WS(w,A,p) for all A, p
% use the value function to update

% guess
WU = 0*JJ;
WU_update = 0*JJ;
dWU2all = 0*JJ';
w  = 1;

JJ = JJ.*(JJ>0);

for x = 1:10
zz = 0    
diffpall = 0*JJ';

for i = 1 : II
%     WUi   = WU(plow_pos(i):end,i)';    
%     diffp = p(plow_pos(i)+1:end) - p(plow_pos(i):end-1); 
%         
%     diffWUx = (WUi(2:end)-WUi(1:end-1))./diffp;
%     dWU2t   = (diffWUx(2:end)-diffWUx(1:end-1))./diffp(1:end-1);
%     dWU2    = [zeros(1,plow_pos(i)) dWU2t dWU2t(end)];
%     dWU2    = dWU2.*(dWU2>-5);

      WUi = WU(:,i)';
      dWU =  (WUi(2:end)-WUi(1:end-1))./(p(2:end)-p(1:end-1));
      ddWU = (dWU(2:end)-dWU(1:end-1))./(p(3:end)-p(2:end-1));
      dWU2 = [0 ddWU ddWU(end)];
      dWU2 = [zeros(1,plow_pos(i)) dWU2(plow_pos(i)+1:end)];
   
      dWU2all(i,:) = dWU2;
%     diffpall(i,plow_pos(i)+1:end) = (p(2:end)-p(1:end-1));
    
%     figure(4)
%     subplot(1,2,1)
%     plot(p,dWU2)
%     subplot(1,2,2)
%     plot(p(plow_pos(i)+1:end),diffWUx)
%     pause      
figure(1)
subplot(1,2,1)
hold off
plot(p(1:end-1),dWU)
hold on
plot(plow(i),dWU(plow_pos(i)),'r.')

subplot(1,2,2)
hold off
plot(p(2:end-1),ddWU)
hold on
plot(plow(i),dWU(plow_pos(i)),'r.')

pause

    
    disp(sprintf('====='))        
    for k = 1 : Np       
    
    if (k > plow_pos(i))    

        %dWU2(plow_pos(i)) = dWU2(plow_pos(i)+1);
%         if dWU2(k)==0 & k ==plow_pos(i)
%         zz = zz+1  ;
%         end
        Term1 = w-b+0.5*sigmap2(k)*dWU2(k);

        % firm-level shock
        Term2 = zeros(1,II);
        Temp2 = zeros(1,4);
        Term3 = zeros(1,II);
        Temp3 = zeros(1,2);
            for j =1:II
                if (j ~= i)
                    pprime = p(k)*(1-Q(i,j))+ HLboth*(1-p(k))*Q(i,j) ;
                    pprimeind = find(p-pprime>0,1);
                        if (pprimeind > 1)
                             wright = (pprime - p(pprimeind-1))/(p(pprimeind)-p(pprimeind-1));
                             new_JJ = (1-wright)*JJ(pprimeind-1,j)+wright*JJ(pprimeind,j);
                             if new_JJ<=0
                                  Term2(j) = -Omega(i,j)*WU(k,i); 
                             elseif new_JJ>0
                                 new_WU = (1-wright)*WU(pprimeind-1,j)+wright*WU(pprimeind,j);
                                 Index2 = 1*(0<new_WU & new_WU<new_JJ) + 2*(new_WU>new_JJ) + 3*(new_WU<=0);
                                 Temp2(1) = new_WU - WU(k,i);
                                 Temp2(2) = new_JJ-WU(k,i);
                                 Temp2(3) = -WU(k,i);
                                 Term2(j) = Omega(i,j)*Temp2(Index2);
                             end                               
                        end
                end
                % search on the job
                Temp3(1) = sum((JJ(k,i) - WU(k,i))*fp1.*dpc/2.*(JJ(:,j)>JJ(k,i))');
                Temp3(2) = sum((JJ(:,j)' - WU(k,i)).*fp1.*dpc/2.*(JJ(:,j)<JJ(k,i) & JJ(:,j)>WU(k,i))');
                Term3(j) = psi*f*vbar(j)*(Temp3(1)+Temp3(2));
            end
            WU_update(k,i) = Term1+ sum(Term2(:)+Term3(:));
    else
        WU_update(k,i) = 0;
    end
    end
end
    figure(2)
    subplot(1,2,1)
    plot(p,WU)
    title('Value function')
    subplot(1,2,2)
    plot(p,dWU2all)
    title('Second derivative')

%     figure(2)
%     subplot(2,1,1)
%     plot(p,dWU2all)
%     subplot(2,1,2)
%     plot(p,diffpall)
%     pause
    
    WU = WU_update/rho;
end


