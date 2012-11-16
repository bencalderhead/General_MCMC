function Result = Evaluate_Prior_Gamma(Paras, Value, Order)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamma prior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Para(1) = parameter 1
% Para(2) = parameter 2

if Order == 0
    % Calculate PDF
    
    if (Value <= 0)
        Result = 0;
    else
        % Calculate probability of value from the prior
        Result = gampdf(Value, Paras(1), Paras(2));
    end
    
    
elseif Order == 1
    % Calculate first derivative of log PDF
    
    if (Value <= 0)
        Result = -1e300;
    else
        Result = (Paras(1)-1)/Value - 1/Paras(2);
    end
    
    
    
elseif Order == 2
    % Calculate second derivative of log PDF
    
    if (Value <= 0)
        Result = -1e300;
    else
        Result = -(Paras(1)-1)/Value^2;
    end
    
    
else
    disp('Error evaluating prior')
end


end

