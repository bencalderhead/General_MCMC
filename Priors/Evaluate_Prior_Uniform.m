function Result = Evaluate_Prior_Uniform(Paras, Value, Order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform prior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Para(1) = Lower bound
% Para(2) = Upper bound

if Order == 0
    
    % Calculate PDF
    
    LowerBound = Paras(1);
    UpperBound = Paras(2);
    
    % Calculate probability of value from the prior
    if (Value > UpperBound) || (Value < LowerBound)
        Result = 0;
    else
        Result = 1/(UpperBound-LowerBound);
    end
    
elseif Order == 1
    % Calculate first derivative of log PDF
    
    LowerBound = Paras(1);
    UpperBound = Paras(2);
    
    if (Value > UpperBound) || (Value < LowerBound)
        Result = -1e300;
    else
        Result = 0;
    end
    
    
    
elseif Order == 2
    % Calculate second derivative of log PDF
    
    LowerBound = Paras(1);
    UpperBound = Paras(2);
    
    if (Value > UpperBound) || (Value < LowerBound)
        Result = -1e300;
    else
        Result = 0;
    end
    
else
    disp('Error evaluating prior')
end


end

