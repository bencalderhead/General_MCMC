function Result = Evaluate_Prior_Normal(Paras, Value, Order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal prior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Para(1) = Mean
% Para(2) = SD

if Order == 0
    
    % Calculate PDF
    
    Mean = Paras(1);
    SD   = Paras(2);
    
    % Calculate probability of value from the prior
    Result = normpdf(Value, Mean, SD);
    
elseif Order == 1
    % Calculate first derivative of log PDF
    
    Mean = Paras(1);
    SD   = Paras(2);
    
    Result = -(Value - Mean)/SD^2;
    
    
    
elseif Order == 2
    % Calculate second derivative of log PDF
    
    Mean = Paras(1);
    SD   = Paras(2);
    
    Result = -1/SD^2;
    
else
    disp('Error evaluating prior')
end


end

