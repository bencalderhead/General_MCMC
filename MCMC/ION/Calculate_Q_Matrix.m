function [Q] = Calculate_Q_Matrix(ModelName, RateParas)

%%% Calculate the Q matrix of rates for given model with specified rate
%%% parameters.

switch ModelName

    case 'Castillo_Katz'
        Q = zeros(3); % 3 states for this model
        
        % Rename for clarity
        K_1   = RateParas(1);
        K_2   = RateParas(2);
        Beta  = RateParas(3);
        Alpha = RateParas(4);
        
        % Construct Q matrix
        Q(1,1) = -K_1;
        Q(1,2) = K_1;
        Q(2,1) = K_2;
        Q(2,2) = -K_2 - Beta;
        Q(2,3) = Beta;
        Q(3,2) = Alpha;
        Q(3,3) = -Alpha;
        
    case 'Five_State'
        % As described in Ch 18 & 20 in Single Channel Recording book
        Q = zeros(5); % 5 states for this model
        
        % Rename for clarity
        K_p1     = RateParas(1);  % k_{+1}
        K_m1     = RateParas(2);  % k_{-1}
        K_p2     = RateParas(3);  % k_{+2}
        K_m2     = RateParas(4);  % k_{-2}
        KStar_p2 = RateParas(5);  % k*_{+2}
        KStar_m2 = RateParas(6);  % k*_{-2}
        Beta_1   = RateParas(7);  % Beta_{1}
        Alpha_1  = RateParas(8);  % Alpha_{1}
        Beta_2   = RateParas(9);  % Beta_{2}
        Alpha_2  = RateParas(10); % Alpha_{2}
        
        
        % Agonist concentration
        X_A = 1e-7; % This should become an extra variable to be passed
        
        % Construct Q matrix - [1 2] are open states, [3 4 5] are closed
        Q(1,1) = -(Alpha_1 + KStar_p2*X_A);
        Q(1,2) = KStar_p2*X_A;
        Q(1,4) = Alpha_1;
        
        Q(2,1) = 2*KStar_m2;
        Q(2,2) = -(Alpha_2 + 2*KStar_m2);
        Q(2,3) = Alpha_2;
        
        Q(3,2) = Beta_2;
        Q(3,3) = -(Beta_2 + 2*K_m2);
        Q(3,4) = 2*K_m2;
        
        Q(4,1) = Beta_1;
        Q(4,3) = K_p2*X_A;
        Q(4,4) = -(Beta_1 + K_p2*X_A +K_m1);
        Q(4,5) = K_m1;
        
        Q(5,4) = 2*K_p1*X_A;
        Q(5,5) = -2*K_p1*X_A;
        
    case 'Five_State_Balanced'
        % As described in Ch 18 & 20 in Single Channel Recording book
        % but with micrscopic reversibility
        Q = zeros(5); % 5 states for this model
        
        % Rename for clarity
        K_p1     = RateParas(1);  % k_{+1}
        K_m1     = RateParas(2);  % k_{-1}
        K_p2     = RateParas(3);  % k_{+2}
        K_m2     = RateParas(4);  % k_{-2}
        Beta_1   = RateParas(5);  % Beta_{1}
        Alpha_1  = RateParas(6);  % Alpha_{1}
        Beta_2   = RateParas(7);  % Beta_{2}
        Alpha_2  = RateParas(8);  % Alpha_{2}
        KStar_p2 = RateParas(9);  % k*_{+2}
        
        % Microsopic reveribility condition
        KStar_m2 = (Alpha_2/Beta_2)*(Beta_1/Alpha_1)*(K_m2/K_p2)*KStar_p2;  % k*_{-2}
        
        
        % Agonist concentration
        X_A = 1e-7; % This should become an extra variable to be passed
        
        % Construct Q matrix - [1 2] are open states, [3 4 5] are closed
        Q(1,1) = -(Alpha_1 + KStar_p2*X_A);
        Q(1,2) = KStar_p2*X_A;
        Q(1,4) = Alpha_1;
        
        Q(2,1) = 2*KStar_m2;
        Q(2,2) = -(Alpha_2 + 2*KStar_m2);
        Q(2,3) = Alpha_2;
        
        Q(3,2) = Beta_2;
        Q(3,3) = -(Beta_2 + 2*K_m2);
        Q(3,4) = 2*K_m2;
        
        Q(4,1) = Beta_1;
        Q(4,3) = K_p2*X_A;
        Q(4,4) = -(Beta_1 + K_p2*X_A +K_m1);
        Q(4,5) = K_m1;
        
        Q(5,4) = 2*K_p1*X_A;
        Q(5,5) = -2*K_p1*X_A;
        
    otherwise
        disp('Error - Invalid Model Name')
end

end
