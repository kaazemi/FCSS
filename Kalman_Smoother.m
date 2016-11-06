function sys_smoothed  = Kalman_Smoother( sys )

% Kalman Smoother for variable
% x_t = A x_{t-1} + w_t
% y_t = C x_t + z_t
% Note: The smoother outputs the augmented estimated vectors if the system
% order is 2, need to downsample by a factor of 2
% For time-varying dynamics simply need to change A -> A(t), C -> C(t) and R -> R(t)
% See Kalman_Filter.m for comments on how to do this

%% Initialization        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = sys.y;                                                              %
[n,T] = size(y);                                                        %    
p = size(sys.C,2);                                                      %
                                                                        %
A = calc_transition_matrix(sys.A,p,sys.Order);                          %
% A = sys.A;                                                              %
                                                                        %
if p ==1                                                                %
Q = reshape(upsample(sys.Q,sys.Order^2),[sys.Order sys.Order, T]);      %    
else                                                                    %
Q = upsample(sys.Q,sys.Order);                                          %
Q = upsample(permute(Q,[2 1 3]),sys.Order);                             %
end                                                                     %
                                                                        %
Xtt = sys.Xtt;                                                          %        
Sigmatt = sys.Sigmatt;                                                  %    
Sigma_st = zeros(size(Sigmatt));                                        %
                                                                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smoothing
    %%%%%%%%%%%%%%%%%%%%% Smoother Initialization %%%%%%%%%%%%%%%%%%%%%%%
    X_smoothed(:,T) = Xtt(:,T);                                         %
    Sigma_smoothed(:,:,T) = Sigmatt(:,:,T);                             %
                                                                        %
    % Calculate Sigma_st(:,:,T) (Required for EM Step)                  %
    t = T;                                                              %
    Sigt = A*Sigmatt(:,:,t)*A'+Q(:,:,t);                                %
        if sys.diagFlagA && sys.diagFlagQ                               %
            St = diag( diag(A).*diag(Sigmatt(:,:,t))./diag(Sigt));      %
        else                                                            %
            St = A*Sigmatt(:,:,t)/Sigt;                                 %
        end                                                             %
    Sigma_st(:,:,t) = St*Sigmatt(:,:,t);                                %
                                                                        %
    %%%%%%%%%%%%%%%%%%%%%%% Backward Iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                  %
    for t=T-1:-1:1                                                                                                %
        Sigt = A*Sigmatt(:,:,t)*A'+Q(:,:,t);                                                                      %  '
        
        % Smoothing Gain (Semi-Optimized)                                                                         %
                St = A*Sigmatt(:,:,t)/Sigt;                                                                       % 
%         switch sys.Order                                                                                          %
%             case 1                                                                                                %
%                 if sys.diagFlagA && sys.diagFlagQ                                                                 %
%                         St = diag( diag(A).*diag(Sigmatt(:,:,t))./diag(Sigt));                                    %
%                 else                                                                                              %
%                         St = A*Sigmatt(:,:,t)/Sigt;                                                               %
%                 end                                                                                               %
%             case 2                                                                                                %
%                         St = [];                                                                                  %
%                 if sys.diagFlagA && sys.diagFlagQ                                                                 %
%                     for i=1:p                                                                                     %
%                         Sti = A(2*i-1:2*i,2*i-1:2*i)*Sigmatt(2*i-1:2*i,2*i-1:2*i,t)/Sigt(2*i-1:2*i,2*i-1:2*i);    %
%                         St = blkdiag(St,Sti);                                                                     %
%                     end                                                                                           %
%                 else                                                                                              %
%                         St = A(:,:,t)*Sigmatt(:,:,t)/Sigt;                                                        %
%                 end                                                                                               %
%         end                                                                                                       %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Update Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                  %
                        X_smoothed(:,t) = Xtt(:,t) + St*(X_smoothed(:,t+1)-A*Xtt(:,t));                           %
                        Sigma_smoothed(:,:,t) = Sigmatt(:,:,t)+St *(Sigma_smoothed(:,:,t+1)-Sigt)*St';            %
                        Sigma_st(:,:,t) = St*Sigma_smoothed(:,:,t);                                               %
    end                                                                                                           %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                  %
sys_smoothed = sys;                                                                                               %
sys_smoothed.X_smoothed = X_smoothed;                                                                             %
sys_smoothed.Sigma_smoothed = Sigma_smoothed;                                                                     %
sys_smoothed.Sigma_st = Sigma_st;                                                                                 %
sys_smoothed.diagFlagA = sys.diagFlagA;                                                                           %
sys_smoothed.diagFlagQ = sys.diagFlagQ;                                                                           %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

