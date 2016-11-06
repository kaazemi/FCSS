function  sys_smoothed = IRLS(sys)

% Iterative Reweighted Least Squares Algorithm for Autoregressive signal deconvolution Model % 
%           x(t) = A(t) x(t-1) + w(t),  w(t) ~ Laplace i.i.d. (lambda)                       %
%           y(t) = C(t) x(t)   + v(t),  v(t) ~ N(0,R(t))                                     %
%                                                                                            %
% Input struct "sys" should include the following fields                                     %
%           sys.y : Observed Calcium Signal, each row is the activity of a single neuron     %
%                                                                                            %
%           sys.lambda: regularization parameter, could choose via cross-validation          %
%                                                                                            %
% Optional: sys.C : The measurement MATRIX C, deafult: Identity (denoising)                  %
%           Should be given for Compressive Calcium Imaging                                  %
%                                                                                            %
%           sys.R : Initial estimate of the observation noise covariance MATRIX              %
%           default: R = sn from FOOPSI (Pnevmatikakis et al.(2016)                          %
%                                                                                            %
%           sys.A : A for a single neuron (scalar or 2 dimensional vector [a b])             %
%           where x(t) = a x(t-1) + b x(t-1) + w(t)                                          %
%           Warning: This will also turn off the EM steps                                    %
%                                                                                            %
%           sys.Order: Autoregressive model order, default: 1                                %
%           default: initialize with aryule(sys.y,sys.Order)                                 %
%                                                                                            %
%           sys.maxNumIters chooses the number of IRLS iterations, default: 5                %
%                                                                                            %
%           sys.resetRFlag if set to true, the estimate of the observation noise covariance  %
%           matrix gets updated, default: false                                              %
%                                                                                            %
%           sys.EMFlag if set to false, EM updates of sys.A will not be                      %
%           performed, default = true                                                        %
%                                                                                            %
%           sys.baseline baseline of calcium measurements, default value is mean of          %
%           anything within 3*std of noise from minimum of traces                            %
%                                                                                            %
%           sys.confidence confidence in detected spikes, default: 90 %                      %
%
%           sys.fs: sampling frequency for spindles; minimum spindle frequency is set to 8Hz %
%           and maximum is set to be 16 Hz                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by:                                                                                %
% Abbas Kazemipour, University of Maryland, College Park,                                    %
%                   Janelia Research Campus         Last Update: November,      1, 2016      %
%                                                                                            %
%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(sys,'Order'); sys.Order = 1; end;                                                %
if ~isfield(sys,'A'); theta = aryule(sys.y',sys.Order);                                      % 
    sys.A = -theta(:,2:end);   end; %sys.EMFlag = false;                                     %
if ~isfield(sys,'confidence'); sys.confidence = 0.90; end                                    %
[n,T] = size(sys.y); % number of measurements per frame (n) and number of time frames (T)    %
if ~isfield(sys,'C'); p = n; sys.C = eye(p);                                                 %
else p = size(sys.C,2);  end                             % State-Space dimension (p);        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The number of iterations is set to be 5 as in most cases this is                     %
% enough for convergence, can easily use convergence criteria (in the main for loop below),  % 
% which is avoided for code simplicity and readability                                       %
maxNumIters = 5;                                                                             %
if isfield(sys,'maxNumIters'); maxNumIters = sys.maxNumIters; end                            %
                                                                                             %
% The parametere epsilon of IRLS should not be chosen too large or too small                 %
% Theoretical results of AR estimation theory suggest very small epsilon lead to             %
% bad estimates of the Model parameteres, whereas large epsilon undesirably                  %
% smoothens the traces, a suggested range is epsilon ~ 1e-6 to 1e-10                         %
epsilon = 1e-15;                                                                             %
                                                                                             %
% lambda can be chosen to be different for different neurons based on the                    %
% estimates of the noise level (modify the code accordingly please)                          %
% lambda is scaled using theoretical guarantees of LASSO                                     %
if isequal(n,p); lambda = sys.lambda * sqrt(T*log(p*T)/n);                                   %
else lambda = sys.lambda * sqrt(log(p*T)/n); end                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate measurement variances                                                            %
if ~isfield(sys,'R'); sys.R = calc_var(sys.y);                                               %
disp('ROI variances calculated from data');    end                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the baseline as the mean of value 3*std of noise above the minimum of traces           %
%                                                                                            %
if isfield(sys,'baseline'); baseline =  sys.baseline;                                        %
else baseline = calc_base(sys);                                                              %
end;                                                                                         %
sys.y = sys.y - repmat(baseline,1,T) ;                                                       %
if ~isfield(sys,'resetRFlag'); sys.resetRFlag = false; end;                                  %
if ~isfield(sys,'EMFlag'); sys.EMFlag = true; end;                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sys.Q, W] = calc_Q(sys.y,sys.A,epsilon,sys.C);                                              %
sys.Q = sys.Q/lambda;                                                                        %
%%%%%%%%%%%%%%%%%%%%% Main IRLS Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iters =1:maxNumIters                                                                     %
    %% Filter Data                                                                           %
    filtered_sys  = Kalman_Filter(sys);                                                      %
    %% Smoothen Data                                                                         %
    sys_smoothed  = Kalman_Smoother(filtered_sys);                                           %
    %% Update State Transition Matrix via EM Algorithm                                       %
    if sys.EMFlag                                                                            %
        sys.A = Expectation_Maximization( sys_smoothed, W );                                 %
        sys_smoothed.A = sys.A;	                                                             %
    end                                                                                      %
    sys = sys_smoothed;                                                                      %
    %% Update the State Covariance Matrix                                                    %
    [sys.Q,W] = calc_Q(sys.X_smoothed(1:sys.Order:end,:),sys.A,epsilon,sys.C);               %
    sys.Q = sys.Q/lambda;                                                                    %
%%%%%%%%%%%%%%%%%%%% Update Estimates of the Noise Variances  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     sys.y = sys.X_smoothed(1:sys.Order:end,:);                                               %
    if sys.resetRFlag && n==p                                                                %
    sys.R = diag( var(sys.y-downsample(sys.X_smoothed,sys.Order),[],2) );                    %
    end                                                                                      %
    %% Replace the initial conditions with the new estimates                                 %
    % Default: x0 = 0; Sig0 = I;                                                             %
    sys.x0 = sys.X_smoothed(:,1);                                                            %
    sys.Sig0 = sys.Sigma_smoothed(:,:,1);                                                    %
end                                                                                          %
                                                                                             %
if sys.Order ==2                                                                             %
sys_smoothed.X_smoothed = downsample(sys_smoothed.X_smoothed,2);                             %
end                                                                                          %
if isequal(n,p)                                                                              %
sys_smoothed.X_smoothed = sys_smoothed.X_smoothed + repmat(baseline,1,T) ;                   %
end                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sys_smoothed.spikes, sys_smoothed.decon, sys_smoothed.decon2] = Find_Spikes(sys_smoothed,sys.confidence);        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

