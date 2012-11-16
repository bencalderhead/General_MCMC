function Result = Sample_Prior_Uniform(Paras)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Uniform prior distribution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Para(1) = Lower bound
        % Para(2) = Upper bound
        
        LowerBound = Paras(1);
        UpperBound = Paras(2);
        
        
        % Produce random value from the prior
        Result = LowerBound + rand*(UpperBound-LowerBound);
        
        
end

