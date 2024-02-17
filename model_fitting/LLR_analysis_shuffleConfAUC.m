% TI control analysis version of mean_dprime_metad.m 

%% now do corrs and plot stuff

%load('MI_LR_1192022_2.mat','dense_I_idx');

for s = 1:subjects
    dolog = 1; % take the log of InhIndx?
    
    AUCconfx = AUCx(s,:); % X or delta model ?
    AUCconfd = AUCd(s,:);
    AUCdec = AUCresp(s,:);
    IDX = dense_I_idx; % *100; % can also check on inh_inxS
    
        % do we take the log?
        if dolog, inhIndx = log10(IDX); else, inhIndx = IDX; end
    
        % plot InhIndex vs AUC for decision vs conf
     % AUC x-model
          f = fit(inhIndx',AUCconfx','poly1');
          fit_valsSim.cx(s,:) = coeffvalues(f);
          fitsSim.cx{s,:} = f;
          
          [r,p] = corr(inhIndx',AUCconfx',"type","Pearson");
          pearsSim.cx(s,:) = [r,p];

          [r,p] = corr(inhIndx',AUCconfx',"type","Spearman");
          spearSim.cx(s,:) = [r,p];
     % AUC d-model
          f = fit(inhIndx',AUCconfd','poly1');
          fit_valsSim.cd(s,:) = coeffvalues(f);
          fitsSim.cd{s,:} = f;
          [r,p] = corr(inhIndx',AUCconfd',"type","Pearson");
          pearsSim.cd(s,:) = [r,p];

          [r,p] = corr(inhIndx',AUCconfd',"type","Spearman");
          spearSim.cd(s,:) = [r,p];

          f = fit(inhIndx',AUCdec','poly1');
          fit_valsSim.d(s,:) = coeffvalues(f);
          fitsSim.d{s,:} = f;
          [r,p] = corr(inhIndx',AUCdec',"type","Pearson");
          pearsSim.d(s,:) = [r,p];

          [r,p] = corr(inhIndx',AUCdec',"type","Spearman");
          spearSim.d(s,:) = [r,p];
end

%% Grab real, biological data from subjects 

cd /Users/shaida/Documents/TI/simulations/complete/
[~, ~, corrs] = corrPlots_sim_shuffleConf();
spear = corrs.spear; pears = corrs.spear; fit_vals = corrs.fit_vals; 

%% now do megan analysis!

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
kde.modelx = ksdensity(data.modelx,xi);% kde.modelx = kde.modelx ./ sum(kde.modelx * .001);
kde.modeld = ksdensity(data.modeld,xi);% kde.modeld = kde.modeld ./ sum(kde.modeld * .001);
kde.human = ksdensity(data.human,xi);% kde.human = kde.human ./ sum(kde.human * .001);

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
figure();
subplot(1,2,1);
Y{:,1} = data.modelx;
Y{:,2} = data.modeld;
Y{:,3} = data.human;

for i = 1:size(Y,2)
[f, u, bb] = ksdensity(Y{i},xi);
 f=f/max(f)*0.3; %normalize
    F(:,i)=f;
    U(:,i)=u;
    MED(:,i)=nanmedian(Y{i});
    MX(:,i)=nanmean(Y{i});
    bw(:,i)=bb;
end

violin(Y,'xlabel',{'x-model','δ-model','humans'},'facecolor',[1 0.75 0;0 0.35 0;.3 .3 .3],'edgecolor','none','bw',bw,'mc','k','medc','r-.')
hold on
yline(0,'k--', 'LineWidth', 1);
ylabel('Fitted slope')
yticks([-0.1 -0.05 0 0.05 0.1])
% legend('Mean', 'Median')

% hold on
% h = hist([data.modelx data.modeld],30); h = h./1000;
% h = bar(h,'grouped');
% set(h(1),'facecolor',[1 0.75 0]);
% set(h(2),'facecolor',[0	0.35 0]);
% title('Models')
% xlabel('Linear fit slope')
% ylabel('Proportion')
% legend('cx','cd')

% subplot(1,3,2);


% subplot(1,3,2);
% hold on
% plot(xi,kde.modelx,'Color',[1 0.75 0],'linewidth',2);
% plot(xi,kde.modeld,'Color',[0 0.35 0],'linewidth',2);
% xline(0,'k--', 'LineWidth', 1);
% xlim([-0.1, 0.1])


% title('Models - KDEs')

subplot(1,2,2);
hist(actualData.LLR_distribution,20)%'facecolor',[.3 .3 .3]);
hold on
xline(0,'k--', 'LineWidth', 1);
xlim([-0.05 0.15])
xlabel('Log likelihood ratio')
ylabel('Frequency')



% subplot(1,3,3);
% hold on
% bins = linspace(-1,1,60);
% h1 = hist(noiseCeiling.LLR_distribution,bins); h1 = h1 ./ 1000;
% h2 = hist(actualData.LLR_distribution,bins); h2 = h2 ./ 20;
% h = bar(bins,[h1; h2]','grouped');
% set(h(1),'facecolor','black');
% set(h(2),'facecolor','red');
% title('Log likelihood ratios');
% xlabel('Log likelihood ratio');
% ylabel('Proportion')
% legend('Noise ceiling - model-based','Empirical data')


% set(gcf,'position',[100,200,1500,500])
makePretty(20)


% Confirm the t-test on the actualData.LLR_distribution
[~,p_LLR,~,stats_LLR] = ttest(actualData.LLR_distribution);

[h_check,p_check,stats_check] = lillietest(actualData.LLR_distribution);

