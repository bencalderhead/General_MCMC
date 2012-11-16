function Result = Sample_Prior_Normal(Paras)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normal prior distribution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Para(1) = Mean
        % Para(2) = SD
        
        Mean = Paras(1);
        SD   = Paras(2);
        
        
        % Produce random value from the prior
        Result = Mean + randn*SD;
      
end

