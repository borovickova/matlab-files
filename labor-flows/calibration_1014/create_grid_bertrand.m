function p = create_grid_bertrand(plow,Np)
                    
    MMM = min(plow);
    pJ = linspace(MMM,1,Np);

    p_add = plow;
    p_step = max(diff(pJ));
    p_new = [];

    for i = 1 : length(p_add)
        [dist,ind] = min(abs(pJ-p_add(i)));
            if (dist < p_step / 100)
                pJ(ind) = p_add(i);
            else
                p_new = [p_new p_add(i)];
            end
    end

    p = unique([pJ p_new]);
end