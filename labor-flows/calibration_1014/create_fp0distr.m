function fp0 = create_fp0distr(pJ,p0,p0_sigma)
    % normal distribution with p0 and p0_sigma
    
    % differences
    %dpc1 = [pJ(2)-pJ(1) 0.5*(pJ(3:end)-pJ(1:end-2)) pJ(end)-pJ(end-1)];

    
    if p0_sigma>0
%         p1 = pdf('Normal',pJ,p0,p0_sigma);
        normalization = normcdf(1,p0,p0_sigma)- normcdf(0,p0,p0_sigma);
%         fp0 = p1/normalization;
        midpoints = [pJ(1) (pJ(1:end-1) + pJ(2:end))/2 pJ(end)];
        quantiles = cdf('Normal',midpoints,p0,p0_sigma);
        fp0 = diff(quantiles)/normalization;
        fp0 = fp0 ./ diff(midpoints);
%         size(fp0)
%         size(midpoints)
%         pause
        
    else
       [a aa] = max(abs(pJ-p0));
       fp0 = zeros(1,length(pJ));
       fp0(aa) = 1; %2/dpc1(aa);
    end

end