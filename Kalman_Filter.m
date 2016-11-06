function filtered_sys = Kalman_Filter(sys)

%% Kalman Filter for time-varying dynamics
% x(t) = A(t) x(t-1) + w(t),  w(t) ~ N(0,Q(t))
% y(t) = C(t) x(t)   + v(t),  v(t) ~ N(0,R(t))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: C should be given in its original form as written above                                      %
% If AR(2) model is being used A should be given in the modified form A = [a1 a2; 1 0];              %            
% If AR(2) model is being used the diagFlagA variable is set to 1 as                                 % 
% default to make parallel processing possible                                                       %
% Q should be the modified form also   Q~ = [Q 0; 0 0];                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If system order is not specified default option is gonna be AR(1)                                  %
if ~isfield(sys,'Order'); sys.Order = 1; disp('Autoregressive(1) model'); end;                       %
                                                                                                     %
% States (Neurons) are assumed independent by default                                                %
diagFlagA = true;  % Independent states                                                              %
diagFlagR = true;  % Independent Measurement Noises                                                  %
diagFlagQ = true;  % Independent Gaussians                                                           %
                                                                                                     %
y = sys.y; [n,T] = size(y);                                                                          %
                                                                                                     %
p = size(sys.C,2);                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If Neurons are independent (i.e. diagFlagA = 1) and working in the                                 % 
% denoising regime (i.e. diagFlagC = 1) huge parallel processing may be achieved                     %
                                                                                                     %
diagFlagC = isdiag(sys.C);                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% More Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = sys.C;                                                                                       %
    if sys.Order == 2                                                                                %
    C = upsample(C',sys.Order)';                                                                     %
    end                                                                                              %
                                                                                                     %
    A = calc_transition_matrix(sys.A,p,sys.Order);                                                   %
                                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismatrix(sys.Q); diagFlagQ = isdiag(sys.Q); sys.Q = repmat(sys.Q,1,1,T); end                      %
if p ==1                                                                                             %
Q = reshape(upsample(sys.Q,sys.Order^2),[sys.Order sys.Order, T]);                                   %    
else                                                                                                 %
Q = upsample(sys.Q,sys.Order);                                                                       %
Q = upsample(permute(Q,[2 1 3]),sys.Order);                                                          %
end                                                                                                  %
%%%%%%%%%%%%%  Optional: Make measurement noise time-varying %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ismatrix(sys.R); diagFlagR = isdiag(sys.R); sys.R = repmat(sys.R,1,1,T); end                    %
R = sys.R;                                                                                           %
if size(R,2)<p; error('R should be diagonal'); end                                                   %
                                                                                                     %
%%%%%%%%%%%%%  Initial Estimates of the States Mean (x0) and Covariance (Sig0)  %%%%%%%%%%%%%%%%%%%%%%
% Note: These may or may not get updated in the IRLS after each Smoothing iteration                  %
                                                                                                     %
if ~isfield(sys,'x0'); sys.x0 = zeros(p*sys.Order,1); end;                                           %
x0 = sys.x0;                                                                                         %
                                                                                                     %
if ~isfield(sys,'Sig0'); sys.Sig0 = eye(p*sys.Order)*mean(diag(R)); end;                             %
Sig0 = sys.Sig0;                                                                                     %       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                     %
Xt = zeros(p*sys.Order,T);                                                                           %
Xtt = zeros(p*sys.Order,T);                                                                          %
Sigt = zeros(p*sys.Order,p*sys.Order,T);                                                             %
Sigmatt = zeros(p*sys.Order,p*sys.Order,T);                                                          %
                                                                                                     %
%% Filtering                                                                                         %
       for t=1:T                                                                                     %
            %% Prediction Step                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Since A is block diagonal under a lot of denoising regimes,                                  % 
% the following multiplications could also be implemented in parallel                                %
                                                                                                     %
            if t==1                                                                                  %
            Xt(:,t) = A*x0;                                                                          %
            Sigt(:,:,t) = A*Sig0*A'+Q(:,:,t);                                                        %
            else                                                                                     %
            Xt(:,t)  = A*Xtt(:,t-1);                                                                 %
            Sigt(:,:,t) = A*Sigmatt(:,:,t-1)*A'+Q(:,:,t);                                            %
            end                                                                                      %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %% Update Step
            
            % Kalman Gain (Semi-Optimized for parallel processing)
                        Kt = (Sigt(:,:,t)*C')/(C*Sigt(:,:,t)*C'+R);
                        
%             switch sys.Order
%                 case 1
%             if diagFlagC && diagFlagA
%                 Kt  = diag( diag(Sigt(:,:,t)).*diag(C)./(diag(Sigt(:,:,t)).*diag(C).^2 + diag(R)) );
%             else
%                 Kt = (Sigt(:,:,t)*C')/(C*Sigt(:,:,t)*C'+R(:,:,t));
%             end
%                 case 2
%                     Kt =[];
%             if diagFlagC && diagFlagA
%                 for i=1:p
%                 Kti = Sigt(2*i-1:2*i,2*i-1:2*i,t)*C(i,2*i-1:2*i)'/( C(i,2*i-1:2*i)*Sigt(2*i-1:2*i,2*i-1:2*i,t)*C(i,2*i-1:2*i)' + R(i,i) );
%                 Kt = blkdiag(Kt,Kti);
%                 end
%             else
%                 Kt = (Sigt(:,:,t)*C')/(C*Sigt(:,:,t)*C'+R);
%             end
%             
%             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: in the denoising regime Kt is block-diagonal and parallel                              %
% processing could potentially make the following multiplications faster                       %
            Xtt(:,t) = Xt(:,t) +Kt*(y(:,t) - C*Xt(:,t));                                       %
            Sigmatt(:,:,t) = Sigt(:,:,t) - Kt*(C*Sigt(:,:,t)*C'+R)*Kt';                        %
       end                                                                                     %
                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Writing the Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only those outputs required for (possible) parallel processing (Flags)                       %
% and Smoothing are given, prediction estimates can also be given if needed                    %
filtered_sys = sys;                                                                            %
filtered_sys.Xtt = Xtt;                                                                        %
filtered_sys.Sigmatt = Sigmatt;                                                                %
% filtered_sys.Xt = Xt;                                                                        %
% filtered_sys.Sigt = Sigt;                                                                    %
%% Flags                                                                                       %
filtered_sys.diagFlagA = diagFlagA;                                                            %
filtered_sys.diagFlagQ = diagFlagQ;                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end