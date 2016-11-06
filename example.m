clear; close all; clc;
%% generate data
T = 1000; nSpikes = floor(T*0.05);
ar_params = [1.4 -.41];
p = 10;
w = zeros(p,T);
for i = 1:p
w(i,randsample(T,nSpikes)) = 1;
end
y = filter(1, [1 -ar_params],w,[],2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sys.y = y + 0.2*randn(size(y));
sys.lambda = 0.1; %Increase for smoother estimates/decrease for noisier

% sys.maxNumIters = 10;      % Increase for better converge (slower) and vice versa
% sys.EMFlag = false;        % can turn EM off
% sys.baseline = zeros(p,1); % can input baseline manually

%% Run IRLS
tic;
sys_smoothed = IRLS(sys);
toc

%% Illustrate
ROI_num = 2;

ax = subplot(3,1,1); hold on
plot(sys.y(ROI_num,:)/max(sys.y(ROI_num,:)),'k');
stem(w(ROI_num,:),'b','linewidth',2,'marker','none');
title('Ground-truth')

bx = subplot(3,1,2); hold on;
plot(sys_smoothed.X_smoothed(ROI_num,:)/max(sys_smoothed.X_smoothed(ROI_num,:)),'k');
% stem(sys_smoothed.spikes(ROI_num,:)/max(sys_smoothed.spikes(ROI_num,:)),'marker','none');
stem(sys_smoothed.spikes(ROI_num,:)>0,'b','linewidth',2,'marker','none');
axis tight
title('FCSS Estimates')

%% If you have foopsi add on path and compare
cx = subplot(3,1,3); hold on;
options.p = 1;
[c,b,c1,g,sn,sp] = constrained_foopsi(y(ROI_num,:),[],[],[],[],options);
plot(c/max(c),'k');
stem(sp>0,'b','linewidth',2,'marker','none');
axis tight
title('Foopsi Estimates')

linkaxes([ax bx cx],'x');
r1 = 0.4; r2 = 0.5;
set(gcf,'units','normalized','outerposition',[0 0 r1 r2],'defaulttextinterpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


