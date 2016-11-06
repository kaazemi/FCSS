function A = Expectation_Maximization(sys_smoothed,W)

% Expectation Maximization Algorithm for estimating the state-transition
% matrix A from the output of a Kalman Smoother
%%%%%%%%%%%%%%%%%%% Model %%%%%%%%%%%%%%%%%%%%%
% x_t = A x_{t-1} + w_t
% y_t = C x_t + z_t
%%%%%%%%%%% Optimization Problem  %%%%%%%%%%%%%
% minimize  Sum E(||x_t - A x_{t-1} ||_1 
%    A      t>1                         
% The problem is solved by IRLS method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Order = sys_smoothed.Order;
X = sys_smoothed.X_smoothed;
Sigma = sys_smoothed.Sigma_smoothed;
Sigma_st = sys_smoothed.Sigma_st;
[p,T] = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma(Sigma<0) = eps;
% Sigma_st(Sigma_st<0) = eps;
% Sigma_st(:,:,end) = Sigma_st(:,:,end) + (Sigma_st(:,:,end-1)>0)*eps;
switch Order
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        mode = 'EM';
        switch mode
            case 'EM'
        Weights = reshape(W(W>0),p,T);        
        X0 = mult_shift(X,X,0); % X_t X_t
        
        X1 = mult_shift(X,X,1); % X_t X_{t-1}

        S = reshape(Sigma,p^2,T);
        
        S0 = S(1:p+1:end,:); % Sigmat-1t-1  
       
        S = reshape(Sigma_st,p^2,T);
        
        S1 = S(1:p+1:end,:); % Sigmatt-1

        num = sum(mult_shift(Weights,X0+S0,1),2);
        denom = sum(mult_shift(Weights,X1+S1,0),2);
        
        A = num./denom;
    
        A = min(A,.999);
    
            case 'yule'
            A = aryule(X(1:2:end,:)',Order);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
    if ~isfield(sys_smoothed,'mode'); mode = 'EM';
    else mode = sys_smoothed.mode;
    end
    
    switch mode
     
        case 'EM'
            %%%%%% Calculate the terms in the equation for EM
        Weights = reshape(W(W>0),p/2,T);
        X = X(1:2:end,:);
        X0 = mult_shift(X,X,0); % X_t X_t
        X1 = mult_shift(X,X,1); % X_t X_{t-1}
        X2 = mult_shift(X,X,2); % X_t X_{t-2}
        
        
        S = reshape(Sigma,p^2,T);
        
        S0 = S(1:2*p+2:end,:); % Sigmat-1t-1  
        
        S1 = S(2:2*p+2:end,:); % Sigmatt-1
        
        S = reshape(Sigma_st,p^2,T);
        
        S2 = S(p+2:2*p+2:end,:); % Sigmatt-2
        
        
        coeffs = zeros(p/2,5);
        coeffs(:,1) = sum(mult_shift(Weights,X0+S0,1),2);
        coeffs(:,2) = sum(mult_shift(Weights,X0+S0,2),2);
        coeffs(:,3) = 2*sum(mult_shift(Weights,X1+S1,1),2);
        coeffs(:,4) = -2*sum(mult_shift(Weights,X1+S1,0),2);
        coeffs(:,5) = -2*sum(mult_shift(Weights,X2+S2,0),2);

        for i=1:p/2
        [A(i,1),A(i,2)] = calc_Calcium_IP(coeffs(i,:));
        end
        
        case 'yule'
            Weights = reshape(W(W>0),p/2,T);
            A = aryule((Weights.*X(1:2:end,:))',Order);

        case 'burg'
            A = arburg((Weights.*X(1:2:end,:))',Order);

            
            
        case 'spindle'
            fs = sys_smoothed.fs;
            min_freq = 12;
            max_freq = 14;
            %%%%%% Calculate the terms in the equation for EM
            %%%%%% x_{t-1} * x_{t-1}
        
        Weights = reshape(W(W>0),p/2,T);
        X = X(1:2:end,:);
        X0 = mult_shift(X,X,0); % X_t X_t
        X1 = mult_shift(X,X,1); % X_t X_{t-1}
        X2 = mult_shift(X,X,2); % X_t X_{t-2}
        
        
        S = reshape(Sigma,p^2,T);
        
        S0 = S(1:2*p+2:end,:); % Sigmat-1t-1  
        
        S1 = S(2:2*p+2:end,:); % Sigmatt-1
        
        S = reshape(Sigma_st,p^2,T);
        
        S2 = S(p+2:2*p+2:end,:); % Sigmatt-2
        
        coeffs = zeros(p/2,5);
        coeffs(:,1) = sum(mult_shift(Weights,X0+S0,1),2);
        coeffs(:,2) = sum(mult_shift(Weights,X0+S0,2),2);
        coeffs(:,3) = 2*sum(mult_shift(Weights,X1+S1,1),2);
        coeffs(:,4) = -2*sum(mult_shift(Weights,X1+S1,0),2);
        coeffs(:,5) = -2*sum(mult_shift(Weights,X2+S2,0),2);
        
        A = zeros(p/2,2);
        for i=1:p/2
        [A(i,1),A(i,2)] = calc_spindle_IP(coeffs(i,:),fs,min_freq,max_freq);
        end
        
    end
    
end


end


function Z = mult_shift(W,X,s)
M = zeros(size(W));
M(:,1:end-s) = W(:,s+1:end);
Z = M.*X;
end
