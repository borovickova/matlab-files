function [random_ini_belief positions] = random_draw_fp0(p_grid,fp0,Nrand)
    
      Npgrid = length(p_grid);
      dpc = [p_grid(2)-p_grid(1) (p_grid(3:end)-p_grid(1:end-2))/2 p_grid(end)-p_grid(end-1)];
    
      Fp0 = cumsum(fp0.*dpc);
      rand_numbers = rand(Nrand,1);
      rand_num_index = sum(repmat(rand_numbers',Npgrid,1)>repmat(Fp0',1,Nrand));
      
      % IF p0 IS BELOW THE LOWEST plow
      % THEN rand_num_index == 0; replace it with 1
      % BUT THE RANDOM_INI_BELIEF IS 0
      random_ini_belief(rand_num_index==0) = -1;
      rand_num_index(rand_num_index==0) = 1;
      random_ini_belief(rand_num_index>0) = p_grid(rand_num_index(rand_num_index>0));
      positions = rand_num_index;

end