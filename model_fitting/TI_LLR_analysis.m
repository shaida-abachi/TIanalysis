% This script allows for comparison between real data and simulated data of
% the Tuned Inhibition model.
%
% Simulated data is the output of fitting the stim to the mean behavior of
% all real subjects, and the tau (post-decisional time) to the mean meta-d'
% of all real subjects 

%% 1. Creating stim (S) variable with fitted S to d' FOR LR RUNS!!!

[lowPE, highPE] = fit_dprime([3 6]); % fit stim to d' of all subjects, 
                                      % with PE:NE ratios of 3 and 6

% [PE6:lowC_lowCoh PE6:lowC_hiCoh PE3:hiC_lowCoh PE3:hiC_hiCoh] % current shape
% hi coh low conflict, PE:NE 6:1 % need to change to this shape
% hi coh hi conflict, PE:NE 3:1
% low coh low conflict, PE:NE 6:1
% low coh hi conflict, PE:NE 3:1

% to change shape from 1x4 to 4x2
stim = [];
stim(1,:) = [highPE(2,2), highPE(2,2)/6]; % 1_6
stim(2,:) = [lowPE(4,2), lowPE(4,2)/3]; % 1_3
stim(3,:) = [highPE(1,2), highPE(1,2)/6]; % 1_6
stim(4,:) = [lowPE(3,2), lowPE(3,2)/3]; % 1_3

all_stimLR = stim;

%% 2. stimulating LR runs with fitted S to d' and tau to meta-d'

subjects = 100;

parfor s = 1:subjects
   [boldLR{s}, perf{s}, correctLR{s}, saveWLR{s}] = simBOLD_iters_tau(all_stimLR, 0.2, 100);
 
    % prepare for AUC
    [bin_levels{s}, activity{s}, response{s}, answer_LR{s}] = LR_sim_output(perf{s}, boldLR{s}, correctLR{s});
end

% calculate AUC
for s = 1:subjects
    for i = 1:20
        sub_activity = activity{s};
        [~,~,~,AUCx(s,i)] = perfcurve(bin_levels{s}.X, sub_activity(:,i),1);
        [~,~,~,AUCd(s,i)] = perfcurve(bin_levels{s}.D, sub_activity(:,i),1);
        [~,~,~,AUCresp(s,i)] = perfcurve(response{s}, sub_activity(:,i),0);
    end
end 

%% 2.1. Plot LR behavioral performance!

makePrettySize = 15; %  size of text in figures
plotLRbehavior(subjects, response, answer_LR, bin_levels, 1, makePrettySize);

%% 3. simulating MI runs with new adjusted stim and tau 

realLRstim = 0.8 * ([0.008571, 0.001429; 0.007500, 0.002500; ...
                0.006857, 0.001143; 0.006000, 0.002000]); % mult by 0.8 because real stim has 0.8 coherence weight 

newStim_mult = all_stimLR ./ realLRstim; % fitted stim/old stim LR

realMIstim = 0.8 * [0.01, 0; 0.007, 0.003; ... 
             0.007, 0; 0.0042, 0.0028]; % mult by 0.8 because real stim has 0.8 coherence weight 

stimMI = newStim_mult .* realMIstim; % multiplicative factor * old stim MI 

% Now, simualte activity in response to MI stim
[boldMI, ~, correctMI, saveWMI] = simBOLD_itersMI(stimMI, 0.2, 100);

% mean simulated BOLD across TRIALS!
lowCd = mean(boldMI(1,:,:),3); % [.01, 0.0] hiD, lowC
hiCd = mean(boldMI(2,:,:),3); % [.006, .004] highC dense
lowCs = mean(boldMI(3,:,:),3); % [.006, 0.0] lowC sparse
hiCs = mean(boldMI(4,:,:),3); % [.0036, .0024] highC sparse

dense_I_idx = (lowCd - hiCd) ./ (lowCd + hiCd); %SPLA changed denominator to sum 05262023
% sparse_I_idx = (lowCs - hiCs) ./ (lowCs + hiCs); % not used, SPLA changed denominator 05262023
% sparse_I_idx = sparse_I_idx * 100;

dense_I_idx = dense_I_idx * 100;
if any(dense_I_idx < 0) % move distribution to above 0
    dense_I_idx = dense_I_idx - min(dense_I_idx) + .0001; % add .0001 so lowest doesn't end up 0, taking log later
end

%% 4. Calculate linear fits corrs for simulated data 

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

%% 5. Grab real, biological data from subjects 

shuffleConf = 0; % shuffle confidence labels to test null model?

cd /Users/shaida/Documents/TI/simulations/complete/
if shuffleConf
    [~, ~, corrs] = corrPlots_sim_shuffleConf();
else
    [~, ~, corrs] = corrPlots_sim();
end

spear = corrs.spear; pears = corrs.spear; fit_vals = corrs.fit_vals;

%% 6. now do megan analysis!

useControl = 0; % do confidence - decisison as a control?

if useControl
    data.human = fit_vals.c(:,1) - fit_vals.d(:,1);
    data.modelx = fit_valsSim.cx(:,1) - fit_valsSim.d(:,1);
    data.modeld = fit_valsSim.cd(:,1) - fit_valsSim.d(:,1);
else
    data.human = fit_vals.c(:,1);
    data.modelx = fit_valsSim.cx(:,1);
    data.modeld = fit_valsSim.cd(:,1);
end


% To visualize data on same scale 
minAll = min([fit_vals.c(:,1); fit_valsSim.cx(:,1); fit_valsSim.cd(:,1)], [], 'all');
maxAll = max([fit_vals.c(:,1); fit_valsSim.cx(:,1); fit_valsSim.cd(:,1)], [], 'all');
% Plot distributions of real, x-model, and delta-model data
figure(); 
subplot(3,1,1); hist(fit_vals.c(:,1)); title('real'); xlim([minAll maxAll]); 
subplot(3,1,2); hist(fit_valsSim.cx(:,1)); title('sim x'); xlim([minAll maxAll]); 
subplot(3,1,3); hist(fit_valsSim.cd(:,1)); title('sim d'); xlim([minAll maxAll]); 


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


% Violin plot  of x-model, delta-model, and real data 
violin(Y,'xlabel',{'x-model','δ-model','humans'},'facecolor',[1 0.75 0;0 0.35 0;.3 .3 .3],'edgecolor','none','bw',bw,'mc','k','medc','r-.')
hold on
yline(0,'k--', 'LineWidth', 1);
ylabel('Fitted slope')
yticks([-0.1 -0.05 0 0.05 0.1])


subplot(1,2,2);
hist(actualData.LLR_distribution,20)%'facecolor',[.3 .3 .3]);
hold on
xline(0,'k--', 'LineWidth', 1);
xlim([-0.5 0.5])
xlabel('Log likelihood ratio')
ylabel('Frequency')

makePretty(makePrettySize)


% Confirm the t-test on the actualData.LLR_distribution
[~,p_LLR,~,stats_LLR] = ttest(actualData.LLR_distribution);

[h_check,p_check,stats_check] = lillietest(actualData.LLR_distribution);
