clear all; close all; clc;

load formegan.mat

% fit_vals.c
% fit_valsSim.cx
% fit_valsSim.cd

useControl = 0;

if useControl
    data.human = fit_vals.c(:,1) - fit_vals.d(:,1);
    data.modelx = fit_valsSim.cx(:,1) - fit_valsSim.d(:,1);
    data.modeld = fit_valsSim.cd(:,1) - fit_valsSim.d(:,1);
else
    data.human = fit_vals.c(:,1);
    data.modelx = fit_valsSim.cx(:,1);
    data.modeld = fit_valsSim.cd(:,1);
end


% Check AUC of model
[~,~,~,AUC_models] = perfcurve([ones(length(data.modelx),1);zeros(length(data.modeld),1)],[data.modelx;data.modeld],0);
% AUC of model is better when we don't use the control, so we will not use
% the control (0.55 vs 0.58)

% Create kdes
xi = -0.25:.001:0.25;
kde.modelx = ksdensity(data.modelx,xi);
kde.modeld = ksdensity(data.modeld,xi);
kde.human = ksdensity(data.human,xi);

% Maximum likelihood ratio we can expect is if we got things PERFECTLY SPOT
% ON, i.e. our human data is the same as the yellow model data
% Find the bins at which we will look up the KDE in this case
for i = 1:length(data.modelx)
    [~,whichBin.modelx(i)] = min(abs(xi - data.modelx(i)));
end

% Now find which bins we will look up for the actual human
for i = 1:length(data.human)
    [~,whichBin.human(i)] = min(abs(xi - data.human(i)));
end

% Look up bins for 'perfect bang on' model and for human
noiseCeiling.modelx = kde.modelx(whichBin.modelx);
noiseCeiling.modeld = kde.modeld(whichBin.modelx);

actualData.modelx = kde.modelx(whichBin.human);
actualData.modeld = kde.modeld(whichBin.human);

% Log likelihood ratios
noiseCeiling.LLR_distribution = log(noiseCeiling.modelx ./ noiseCeiling.modeld);
actualData.LLR_distribution = log(actualData.modelx ./ actualData.modeld);

noiseCeiling.LLR_mean = mean(noiseCeiling.LLR_distribution);
actualData.LLR_mean = mean(actualData.LLR_distribution);


% Do some traditional stats
[~,p,ci,stats] = ttest(data.human);

% Figures
figure(1);
subplot(1,3,1);
hold on
h = hist([data.modelx data.modeld],30); h = h./1000;
h = bar(h,'grouped');
set(h(1),'facecolor','y');
set(h(2),'facecolor','g');
title('Models')
xlabel('Linear fit slope')
ylabel('Proportion')
legend('cx','cd')

subplot(1,3,2);
hold on
plot(xi,kde.modelx,'y','linewidth',2);
plot(xi,kde.modeld,'g','linewidth',2);
title('Models - KDEs')

subplot(1,3,3);
hold on
bins = linspace(-1,1,60);
h1 = hist(noiseCeiling.LLR_distribution,bins); h1 = h1 ./ 1000;
h2 = hist(actualData.LLR_distribution,bins); h2 = h2 ./ 20;
h = bar(bins,[h1; h2]','grouped');
set(h(1),'facecolor','black');
set(h(2),'facecolor','red');
title('Log likelihood ratios');
xlabel('Log likelihood ratio');
ylabel('Proportion')
legend('Noise ceiling - model-based','Empirical data')


set(gcf,'position',[100,200,1500,500])
makePretty(15)


% Confirm the t-test on the actualData.LLR_distribution
[~,p_LLR,~,stats_LLR] = ttest(actualData.LLR_distribution);

