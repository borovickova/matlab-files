function surplus = calculate_surplus(p,JJ,current_belief,current_firm_shock)

    II = size(JJ,2);
    Nx = length(current_belief);
    surplus = zeros(1,Nx);
    
    for i=1:II
        tempi = current_belief(current_firm_shock==i);
        surplus(current_firm_shock==i) = interp1(p,JJ(:,i)',tempi);
    end
end