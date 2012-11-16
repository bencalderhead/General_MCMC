function Result = Sample_Prior_Gamma(Paras)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gamma prior distribution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Para(1) = parameter 1
        % Para(2) = parameter 2
        
        
        % Produce random value from the prior
        Result = gamrnd(Paras(1), Paras(2));
        

end

