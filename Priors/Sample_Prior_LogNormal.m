function Result = Sample_Prior_LogNormal(Paras)

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Log-Normal prior distribution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Para(1) = Mean
        % Para(2) = SD
        
        Mean = Paras(1);
        SD   = Paras(2);
        
        % Produce random value from the prior
        Result = lognrnd(Mean,SD);


end

