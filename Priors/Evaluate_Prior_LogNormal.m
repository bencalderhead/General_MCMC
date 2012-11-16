function Result = Evaluate_Prior_LogNormal(Paras, Value, Order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log-Normal prior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Para(1) = Mean
% Para(2) = SD

if Order == 0
    % Calculate PDF
    
    Mean = Paras(1);
    SD   = Paras(2);
    
    % Calculate probability of value from the prior
    Result = lognpdf(Value, Mean, SD);
    
elseif Order == 1
    % Calculate second derivative of log PDF
    
    Mean = Paras(1);
    SD   = Paras(2);
    
    Result = -(1/Value) - ((log(Value)-Mean)/(SD^2))*(1/Value);
    
    
    
    
elseif Order == 2
    % Calculate second derivative of log PDF
    
    Mean = Paras(1);
    SD   = Paras(2);
    
    Result = (1/Value^2) + ((log(Value)-Mean)/(SD^2))*(1/Value^2) - (1/SD^2)*(1/Value^2);
    
    
else
    disp('Error evaluating prior')
end



end

