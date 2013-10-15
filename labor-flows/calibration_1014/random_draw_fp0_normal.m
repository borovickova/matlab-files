function [rand_numbers] = random_draw_fp0_normal(Nrand,p0,p0_sigma)

Nrandmax = normcdf(1,p0,p0_sigma);
Nrandmin = normcdf(0,p0,p0_sigma);
rand_numbers = rand(Nrand,1)*(Nrandmax-Nrandmin) + Nrandmin;
rand_numbers = norminv(rand_numbers,p0,p0_sigma);
     
end