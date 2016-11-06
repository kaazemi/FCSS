function [spikes,decon,decon2]  = Find_Spikes( sys, confidence_level )

% Inputs: structure sys: Output of the IRLS method                                     %
% confidence_level: Desired confidence level, a scalar in (0,1), default: 0.9          %
%%%%%%%%%%%%%%% Set default confidence level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2; confidence_level = 0.9; end                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate of how big the error will get relatively due to the differencing            %      
% operator, always set it to 1.                                                        %
alpha = 1;                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p, T] = size(sys.X_smoothed);                                                         %
n = size(sys.y,1);                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deconvolved Signal                                                                   %
decon = zeros(p,T);                                                                    %
for i = 1:p                                                                            %
    decon(i,:) = filter([1 -sys.A(i,:)],1,sys.X_smoothed(i,:),[],2);                   %
    decon2(i,:) = filter([1 -1.999 1],1,sys.X_smoothed(i,:),[],2);                     %
end                                                                                    %
decon(:,1)   = 0; decon(decon<0) = 0;                                                    %
decon2(:,1:2)= 0; decon2(decon2<0) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch n                                                                               %
    case p                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The k-th element of confidence is the confidence bound for the k-th Neuron           %
% See [Van De Geer 2014] for the derivation                                            %
confidence = 1/2*norminv(1-(1-confidence_level)/2)*sqrt(diag(sys.R));                  %
confidence = repmat(confidence,1,T);                                                   %
spikes = decon;                                                                        %
spikes(spikes<confidence) = 0;                                                         %
spikes(:,1:end-2) = spikes(:,3:end); spikes(:,end) = 0; 
spikes(find_dec(sys.X_smoothed)) = 0;
% spikes(find_dec(decon)) = 0;
decon2(spikes == 0) = 0;
for i=1:p
spikes(i,:) = merge_spikes(spikes(i,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    otherwise                                                                          %
        disp('Confidence bounds for Compressive Regime Will be added soon')            %
% Augmented_meas_mtx = toeplitz(ones(T,1),[1;zeros(T-1,1)]);                           
% Augmented_Cov = Augmented_meas_mtx'*Augmented_meas_mtx/T;
% Bias =  (Augmented_Cov\(Augmented_meas_mtx*(sys.y-X_smoothed)')/T)';                 %
% X_smoothed_unbiased = sys.deconv + Bias;                                             %
% X_smoothed_unbiased = X_smoothed + Bias*Augmented_meas_mtx';                         %
% Compensate for the bias caused by regularization                                     %
end                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end


function [X_dec] = find_dec(X)
%%%% Locations where calcium has decreased
X_dec = filter([1 -1],1,X,[],2)<0;
%%%% Locations where calcium has increased
% X_inc = ~X_dec;
% X_dec = X_dec;
% X_inc = X_inc.*X;
end





% function [pks, loc_pks] = findpeaks_mat(X)
% [p,T] = size(X);
% S = vec(X'); 
% [vec_pks, loc_pks] = findpeaks(S);
% pks = zeros(p*T,1);
% pks(loc_pks) = S(loc_pks); pks(1:T:end) = 0;
% pks = reshape(pks,T,p)';
% 
% end

% function X_th = threshold_passed(X,th)
% thresholds = repmat(th,1,size(X,2));
% X_th_flag = X>thresholds;
% 
% filter([1 -1],1,X_th_flag,[],2) == 1;
% 
% end
% 
% 
