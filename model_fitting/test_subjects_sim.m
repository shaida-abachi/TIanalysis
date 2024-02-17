
% 05152023: this script runs sim_MI_LR_new but for many subjects. first
% made with the intention of creating histograms of correlation
% coefficients to then compare to real data

subjects = 1;

% parfor s = 1:subjects
%     [boldMI{s}, ~, correctMI{s}, saveWMI{s}] = simBOLD_iters2('MI', 0.2, 100);
%     sub_boldMI = boldMI{s};
%     
%     lowCd = mean(sub_boldMI(1,:,:),3); % [.01, 0.0] hiD, lowC MEANING ACROSS TRIALS!
%     hiCd = mean(sub_boldMI(2,:,:),3); % [.006, .004] highC dense
%     lowCs = mean(sub_boldMI(3,:,:),3); % [.006, 0.0] lowC sparse
%     hiCs = mean(sub_boldMI(4,:,:),3); % [.0036, .0024] highC sparse
%     
%     dense_I_idx(s) = (lowCd - hiCd) ./ lowCd;
%     sparse_I_idx(s) = (lowCs - hiCs) ./ lowCs;
% end

for s = 2
%   [PE6:lowC_lowCoh PE6:lowC_hiCoh PE3:hiC_lowCoh PE3:hiC_hiCoh]
% hi coh low conflict, PE:NE 6:1
 % hi coh hi conflict, PE:NE 3:1
% low coh low conflict, PE:NE 6:1
% low coh hi conflict, PE:NE 3:1
    stim(1,:) = [PE6{s,2}(2), PE6{s,2}(2)/6]; 
    stim(2,:) = [PE3{s,4}(2), PE3{s,4}(2)/3];
    stim(3,:) = [PE6{s,1}(2), PE6{s,1}(2)/6];
    stim(4,:) = [PE3{s,3}(2), PE3{s,3}(2)/3];
%     all_stim{s} = stim; 
end

parfor s = 1:subjects
    tic
%     waitbar(sub/10,h)
%     [boldLR{s}, perf{s}, correctLR{s}, saveWLR{s}] = simBOLD_iters2('LR', 0.2, 100);
%     [boldLR{s}, perf{s}, correctLR{s}, saveWLR{s}] = simBOLD_iters2(all_stim{s}, 0.2, 100);
      [boldLR{s}, perf{s}, correctLR{s}, saveWLR{s}] = simBOLD_iters2(stim, 0.2, 1000);
    
    % prepare for AUC
    [bin_levels{s}, activity{s}, response{s}, answer_LR{s}] = LR_sim_output(perf{s}, boldLR{s}, correctLR{s});
    time(s) = toc;
end

% calculate AUC
for s = 1:subjects
%     for i = 1:size(boldLR,2)
    for i = 1:20
        sub_activity = activity{s};
        [~,~,~,AUCx(s,i)] = perfcurve(bin_levels{s}.X, sub_activity(:,i),1);
        [~,~,~,AUCd(s,i)] = perfcurve(bin_levels{s}.D, sub_activity(:,i),1);
        [~,~,~,AUCresp(s,i)] = perfcurve(response{s}, sub_activity(:,i),0);
    end
end 

%AUCx, AUCd, AUCresp is subjects x inh interneuron 

% DON'T FORGEGT TO SAVE IT !!!!!!
%save(['' subjects '_subjects_sim.mat'])

%% getting behavioral performance stats and plotting them
% 
% stim = [ 0.008571, 0.001429; ... % hi coh low conflict, PE:NE 6:1
%          0.007500, 0.002500; ...   % hi coh hi conflict, PE:NE 3:1
%          0.006857, 0.001143; ... % low coh low conflict, PE:NE 6:1
%          0.006000, 0.002000 ]; % low coh hi conflict, PE:NE 3:1
 
% variables are set up so out the 400, each increment of 100 is a
 % different stim, so need to separate them by stim:
% stim = [1 100; 101 200; 201 300; 301 400];
stim = [1 1000; 1001 2000; 2001 3000; 3001 4000];

ds = nan(4,subjects);
mcs_x = ds;
mcs_d = ds;
pcs = ds;

for i_sub = 1:subjects
    for i_stim = 1:4
      
       resp = response{i_sub}; % grab data from current simulated subject
       resp = resp(stim(i_stim,1):stim(i_stim,2)); % now grab the correct stim combination (1:4)
       
       correct_ans = answer_LR{i_sub};
       correct_ans = correct_ans(stim(i_stim,1):stim(i_stim,2));
      
       conf_x = bin_levels{i_sub}.X; 
       conf_x = conf_x(stim(i_stim,1):stim(i_stim,2)) + 1;
      
       conf_d = bin_levels{i_sub}.D;
       conf_d = conf_d(stim(i_stim,1):stim(i_stim,2)) + 1;

        % calculate percent correct and mean confidence 
        pc(i_stim) = mean(resp == correct_ans);
        mc_x(i_stim) = mean(conf_x);
        mc_d(i_stim) = mean(conf_d);
        
        % calculate d'
        hr = (sum(correct_ans == 1 & resp == 1) +.5) / (sum(correct_ans == 1)+1);
        far = (sum(correct_ans == 0 & resp == 1) +.5) / (sum(correct_ans == 0)+1);

        d(i_stim)  = norminv(hr) - norminv(far);
        
        % coherence
%             coh(i_stim) = unique(data.stim.fractionCoherent(f));
    end
    ds(:,i_sub)   = d;
    mcs_x(:,i_sub)  = mc_x;
    mcs_d(:,i_sub)  = mc_d;
    pcs(:,i_sub)  = pc;
end

% compute mean and SEM
d_m   = mean(ds,2);
mc_m_x  = mean(mcs_x,2);
mc_m_d  = mean(mcs_d,2);
pc_m  = mean(pcs,2);

d_sem   = std(ds,[],2) / sqrt(subjects);
mc_sem_x  = std(mcs_x,[],2) / sqrt(subjects);
mc_sem_d  = std(mcs_d,[],2) / sqrt(subjects);
pc_sem  = std(pcs,[],2) / sqrt(subjects);

%reshaping things to be like real behavioral data
d_m2(:,1) = d_m(3:4); d_m2(:,2) = d_m(1:2);
mc_m_x2(:,1) = mc_m_x(3:4); mc_m_x2(:,2) = mc_m_x(1:2);
mc_m_d2(:,1) = mc_m_d(3:4); mc_m_d2(:,2) = mc_m_d(1:2);
pc_m2(:,1) = pc_m(3:4); pc_m2(:,2) = pc_m(1:2);

d_sem2(:,1) = d_sem(3:4); d_sem2(:,2) = d_sem(1:2);
mc_sem_x2(:,1) = mc_sem_x(3:4); mc_sem_x2(:,2) = mc_sem_x(1:2);
mc_sem_d2(:,1) = mc_sem_d(3:4); mc_sem_d2(:,2) = mc_sem_d(1:2);
pc_sem2(:,1) = pc_sem(3:4); pc_sem2(:,2) = pc_sem(1:2);

% bar plot d' and mean conf
figure; 

% subplot(1,3,1);
% bar_error(coh_m, coh_sem);
% set(gca,'XTick',[1 2])
% set(gca,'XTickLabel',{'high PE:NE','low PE:NE'})
% ylabel('motion coherence')


subplot(1,3,1);
bar_error(d_m2, d_sem2);
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'high PE:NE','low PE:NE'})
ylabel('d''')
title(['data for n = ' num2str(subjects)])

subplot(1,3,2);
bar_error(mc_m_x2, mc_sem_x2);
yy = get(gca, 'ylim');
set(gca, 'ylim', [1 yy(2)]);
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'high PE:NE','low PE:NE'})
legend('low coherence','high coherence')
ylabel('mean conf x model')

subplot(1,3,3);
bar_error(mc_m_d2, mc_sem_d2);
yy = get(gca, 'ylim');
set(gca, 'ylim', [1 yy(2)]);
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'high PE:NE','low PE:NE'})
legend('low coherence','high coherence')
ylabel('mean conf delta model')

% plot conf_x vs d'
figure; hold on;
plot(d_m2(1,:), mc_m_x2(1,:), 'go-', 'LineWidth', 2);
plot(d_m2(2,:), mc_m_x2(2,:), 'ro-', 'LineWidth', 2);
hv_errorbars(d_m2(1,:), mc_m_x2(1,:), d_sem2(1,:), mc_sem_x2(1,:), 'g-', .5);
hv_errorbars(d_m2(2,:), mc_m_x2(2,:), d_sem2(2,:), mc_sem_x2(2,:), 'r-', .5);
% xlim([0.8, 2.8]);
xlabel('d''')
ylabel('mean conf x model')
title(['data for n = ' num2str(subjects)])

% legend(['low conflict (PE/NE = ' num2str(PE_NE_low) ')'], ['high conflict (PE/NE = ' num2str(PE_NE_high) ')'],'Location','NorthWest');
legend('low PE', 'high PE')

% plot conf_d vs d'
figure; hold on;
plot(d_m2(1,:), mc_m_d2(1,:), 'go-', 'LineWidth', 2);
plot(d_m2(2,:), mc_m_d2(2,:), 'ro-', 'LineWidth', 2);
hv_errorbars(d_m2(1,:), mc_m_d2(1,:), d_sem2(1,:), mc_sem_d2(1,:), 'g-', .5);
hv_errorbars(d_m2(2,:), mc_m_d2(2,:), d_sem2(2,:), mc_sem_d2(2,:), 'r-', .5);
% xlim([0.8, 2.8]);
xlabel('d''')
ylabel('mean conf d model')
title(['data for n = ' num2str(subjects)])

% legend(['low conflict (PE/NE = ' num2str(PE_NE_low) ')'], ['high conflict (PE/NE = ' num2str(PE_NE_high) ')'],'Location','NorthWest');
legend('low PE', 'high PE')

%% now do corrs and plot stuff

load('MI_LR_1192022_2.mat','dense_I_idx');

for s = 1:subjects
    dolog = 1; % take the log of InhIndx?
    
    AUCconfx = AUCx(s,:); % X or delta model ?
    AUCconfd = AUCd(s,:);
    AUCdec = AUCresp(s,:);
    IDX = dense_I_idx*100; % can also check on inh_inxS
    
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


[~, ~, corrs] = corrPlots();
spear = corrs.spear; pears = corrs.spear; fit_vals = corrs.fit_vals;

%% plot corrs 
nbins = 25;

figure(99); clf
subplot(3,3,4)
histogram(atanh(spearSim.cx(:,1)),nbins,'FaceColor','y')
hold on
histogram(atanh(spearSim.cd(:,1)),nbins,'FaceColor','g')
title('spearman')
ylabel('simulated stuff')

subplot(3,3,5)
histogram(atanh(pearsSim.cx(:,1)),nbins,'FaceColor','y')
hold on
histogram(atanh(pearsSim.cd(:,1)),nbins,'FaceColor','g')
title('pearson')

subplot(3,3,6)
histogram(fit_valsSim.cx(:,1),nbins,'FaceColor','y')
hold on
histogram(fit_valsSim.cd(:,1),nbins,'FaceColor','g')
title('linear fit')


% %measures of separation - AUC
% labels = [zeros(100,1); ones(100,1)];
% xvals.s = [atanh(spearSim.cx(:,1)); atanh(spearSim.cd(:,1))];
% xvals.p = [atanh(pearsSim.cx(:,1)); atanh(pearsSim.cd(:,1))];
% xvals.lin = [fit_valsSim.cx(:,1); fit_valsSim.cd(:,1)];
% 
% [~,~,~,histAUC.s] = perfcurve(labels,xvals.s,1);
% [~,~,~,histAUC.p] = perfcurve(labels,xvals.p,1);
% [~,~,~,histAUC.lin] = perfcurve(labels,xvals.lin,1);

%% KDE

% need to do below so atanh doesn't become infinity 
spearSim.cx(spearSim.cx == -1) = -0.9999999; spearSim.cd(spearSim.cd == -1) = -0.9999999;
spearSim.cx(spearSim.cx == 1) = 0.9999999; spearSim.cd(spearSim.cd == 1) = 0.9999999;
figure(99);
dx = 0.001;
xi = -.5:dx:.5;

% SPEARMAN
fSpear.cx = ksdensity(atanh(spearSim.cx(:,1)),xi); fSpear.cx = fSpear.cx ./ sum(fSpear.cx);% bad because infinity, changed 1's to 0.99999 above 
fSpear.cd = ksdensity(atanh(spearSim.cd(:,1)),xi); fSpear.cd = fSpear.cd ./ sum(fSpear.cd); % normalizing it 
bio_spear = round(atanh(spear.c(:,1)),3);
fSpear.empirical = ksdensity(bio_spear,xi); fSpear.empirical = fSpear.empirical ./ sum(fSpear.empirical);
for s = 1:20
    [~,binS(s)] = min(abs(bio_spear(s)-xi)); % find the height of where real data and KDE of sim meet
end
% look up the KDE height at each of the biodata points for yellow & green
px_yellow.spear = fSpear.cx(binS); px_green.spear = fSpear.cd(binS);
BF_spear = mean(log(px_yellow.spear ./ px_green.spear));

subplot(3,3,1); hold on % plot KDE distributions and real data 
plot(xi,fSpear.cx,'y','linewidth',2)
plot(xi,fSpear.cd,'g','linewidth',2)
plot(xi,fSpear.empirical,'k','linewidth',2)
title('Spearman')

subplot(3,3,7); hold on
histogram(log(px_yellow.spear ./ px_green.spear),nbins,'FaceColor','k')
ylabel('actual bayes factors')

% PEARSON
fPears.cx = ksdensity(atanh(pearsSim.cx(:,1)),xi); fPears.cx = fPears.cx ./ sum(fPears.cx);
fPears.cd = ksdensity(atanh(pearsSim.cd(:,1)),xi); fPears.cd = fPears.cd ./ sum(fPears.cd);
bio_pears = round(atanh(pears.c(:,1)),3);
fPears.empirical = ksdensity(bio_pears,xi); fPears.empirical = fPears.empirical ./ sum(fPears.empirical);
for s = 1:20
    [~,binP(s)] = min(abs(bio_pears(s)-xi)); % find the height of where real data and KDE of sim meet
end
% look up the KDE height at each of the biodata points for yellow & green
px_yellow.pears = fPears.cx(binP); px_green.pears = fPears.cd(binP);
BF_pears = mean(log(px_yellow.pears ./ px_green.pears));

subplot(3,3,2); hold on % plot KDE distributions and real data 
plot(xi,fPears.cx,'y','linewidth',2)
plot(xi,fPears.cd,'g','linewidth',2)
plot(xi,fPears.empirical,'k','linewidth',2)
title('Pearson')
ylabel('probability')
subplot(3,3,8); hold on
histogram(log(px_yellow.pears ./ px_green.pears),nbins,'FaceColor','k')

% LINEAR FIT 
fFit.cx = ksdensity(fit_valsSim.cx(:,1),xi); fFit.cx = fFit.cx ./ sum(fFit.cx);
fFit.cd = ksdensity(fit_valsSim.cd(:,1),xi); fFit.cd = fFit.cd ./ sum(fFit.cd);
bio_fit = round(fit_vals.c(:,1),3);
fFit.empirical = ksdensity(bio_fit,xi); fFit.empirical = fFit.empirical ./ sum(fFit.empirical);
for s = 1:20
    [~,binF(s)] = min(abs(bio_fit(s)-xi)); % find the height of where real data and KDE of sim meet
end
% look up the KDE height at each of the biodata points for yellow & green
px_yellow.fit = fFit.cx(binF); px_green.fit = fFit.cd(binF);
BF_fit = mean(log(px_yellow.fit ./ px_green.fit));

%test
TESTbio_fitLEFT = round(fit_vals.c(:,1),3) - (mean(round(fit_vals.c(:,1),3)) - sum(fFit.cx .* xi));
TESTbio_fitRIGHT = round(fit_vals.c(:,1),3) - (mean(round(fit_vals.c(:,1),3)) - sum(fFit.cd .* xi));
for s = 1:20
    [~,TESTbinF_LEFT(s)] = min(abs(TESTbio_fitLEFT(s)-xi)); % find the height of where real data and KDE of sim meet
    [~,TESTbinF_RIGHT(s)] = min(abs(TESTbio_fitRIGHT(s)-xi));
end
TESTpx_leftShiftYellow.fit = fFit.cx(TESTbinF_LEFT); TESTpx_leftShiftGreen.fit = fFit.cd(TESTbinF_LEFT);
TESTpx_rightShiftYellow.fit = fFit.cx(TESTbinF_RIGHT); TESTpx_rightShiftGreen.fit = fFit.cd(TESTbinF_RIGHT);

MAX_BF_fit = mean(log(TESTpx_leftShiftYellow.fit ./ TESTpx_leftShiftGreen.fit));
MIN_BF_fit = mean(log(TESTpx_rightShiftYellow.fit ./ TESTpx_rightShiftGreen.fit));
maxtest = BF_fit/MAX_BF_fit;
mintest =  BF_fit/MIN_BF_fit;
%end test

subplot(3,3,3); hold on % plot KDE distributions and real data 
plot(xi,fFit.cx,'y','linewidth',2)
plot(xi,fFit.cd,'g','linewidth',2)
plot(xi,fFit.empirical,'k','linewidth',2)
xlim([-.1,.1])
title('Linear fit')

subplot(3,3,9); hold on % plot histogram of real data
histogram(log(px_yellow.fit ./ px_green.fit),nbins,'FaceColor','k')
