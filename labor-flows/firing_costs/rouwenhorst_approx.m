
function [Omega AA] = rouwenhorst_approx(quar_rhoy,quar_sigmae,II)
%local monthly_rhoy monthly_sigmae pp qq monthly_sigmay Amax Q1 Q2 Q3 Q4 QQ xx diag1x

% quar_rhoy = 0.95;
% quar_sigmae = 0.3;
% II = 5;

monthly_rhoy = quar_rhoy^(1/3);
monthly_sigmae = quar_sigmae/sqrt(1+quar_rhoy^2+quar_rhoy^4);

pp = 0.5*(1+monthly_rhoy);
qq = pp;
monthly_sigmay = (monthly_sigmae^2/(1-monthly_rhoy^2))^0.5;

Amax = sqrt(II-1)*monthly_sigmay;
AA = linspace(-Amax,Amax,II);
%AA = linspace(-0.2,0.4,II);

QQ = [pp 1-pp; 1-qq qq];
for i=3:II
    Q1 = zeros(i,i);
    Q2 = zeros(i,i);
    Q3 = zeros(i,i);
    Q4 = zeros(i,i);
    
    Q1(1:i-1,1:i-1) = QQ;
    Q2(1:i-1,2:i) = QQ;
    Q3(2:i,1:i-1) = QQ;
    Q4(2:i,2:i) = QQ;

    QQ = pp*Q1+(1-pp)*Q2+(1-qq)*Q3+qq*Q4;
    QQ(2:end-1,:)=QQ(2:end-1,:)/2;
end

Omega = logm(QQ);
Omega = round(Omega*1000)/1000;

diag1x = sum(Omega,2)-diag(Omega);
xx = [0:1:II-1]*II + [1:1:II];
Omega(xx) = -diag1x;

end