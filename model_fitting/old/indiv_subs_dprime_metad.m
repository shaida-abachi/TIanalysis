
% trying individual subjects fitted tau with alpha = 0.1

% % 0.3
% fit_3 = [-457.295932877911	98.0204755844046	0.0103072194984300];
% 
% % 0.6
% fit_6 = [-602.813230649464	121.561575089178	0.0176512639212762];

% 0.1 
fit_1 = [-879.941967991654 136.383362056827	-0.0103557840141724];

% [PE6:lowC_lowCoh PE6:lowC_hiCoh PE3:hiC_lowCoh PE3:hiC_hiCoh]
 sub(1,:)= [2.58046643 2.872464169 1.319627678 2.4148281];
 sub(2,:)= [0.6243293747 0.4826132377 0.3744850743 0.1452834641];
 sub(3,:)= [1.043242296 3.157590626 0.505668042 1.376911035];
 sub(4,:)= [1.336925018 1.78575868 0.9065278437 1.49233436];
 sub(5,:)= [3.039179574 4.482805455 1.155937203 3.597714473];
 sub(6,:)= [1.932778776 3.048378033 1.159111084 2.411120151];
 sub(7,:)= [2.569651368 4.443039177 1.727277381 3.560928683];
 sub(8,:)= [1.698602401 2.30495999 0.9812839733 1.854568125];
 sub(9,:)= [1.920058123 3.169084264 0.9880277804 2.446372699];
 sub(10,:)= [0.1342400125 1.843932025 0.1734856764 0.5162146335];
 sub(11,:)= [0.5621005216 1.95091954 0.306910977 1.704888534];
 sub(12,:)= [1.04022067 0.7982484692 0.6249876153 0.530232748];
 sub(13,:)= [0.2903511651 0.5601963875 -0.1225625209 -0.294534835];
 sub(14,:)= [2.63314912 3.04574771 1.565672481 2.154791089];
 sub(15,:)= [2.030232308 2.375094433 1.442148721 2.365025102];
 sub(16,:)= [2.03272232 2.732382615 1.237513809 1.953568025];
 sub(17,:)= [1.577964419 1.766625025 1.038343286 1.700674488];
 sub(18,:)= [2.333901174 2.834820325 1.25560042 2.066470528];
 sub(19,:)= [0.8619007826 1.005180297 0.3350013948 0.4666264737];
 sub(20,:)= [1.480967952 2.877529792 3.136776087 3.476111938];

 all_lowPE_dp = mean(sub(:,1:2),2);
%  all_dp = mean(sub,1);

for s = 1:20
    for stim = 1:4
        p1 = fit_1;
        p1(3) = fit_1(3) -  sub(s,stim);
        roots_sub = roots(p1); 
        PE01(s,stim) = roots_sub(2);
    end
    lowPE_all = fit_1;
    lowPE_all(3) = fit_1(3) - all_lowPE_dp(s);
    roots_lowPE = roots(lowPE_all); 
    all_lowPE01(s) = roots_lowPE(2);
end

%sim.S = 

% for stim = 1:4
%     p3 = fit_1;
%     p3(3) = fit_1(3) -  all_dp(stim);
%     PE3(stim,:) = roots(p3);
% 
% %     p6 = fit_6;
% %     p6(3) = fit_6(3) -  all_dp(stim);
% %     PE6(stim,:) = roots(p6);
% end
% 
% p3 = fit_1;
% p3(3) = fit_1(3) -  mean(sub(:,3:4),'all');
% lowPE3 = roots(p3);
% 
% p6 = fit_6;
% p6(3) = fit_6(3) -  mean(sub(:,3:4),'all');
% lowPE6 = roots(p6);



%%
subjectIDs = {'101','105','110','111','203','204','205','206','207','208'.../
    ,'209','213','215','217','218','220','221','222','223','224'};
% subjectIDs = {'101'};

for sub = 1:length(subjectIDs)
    data = TN_fMRI_getSubjectData(subjectIDs{sub}, 'LR');
    conf_1(sub) = sum(data.conf==1)/sum(data.conf>0);
    conf_2(sub) = sum(data.conf==2)/sum(data.conf>0);
 % meta d'
%     low_PE = data.stim.PE_NE_ratio == 3;
    low_conflict = data.stim.conflictID == 1;
%     [nR_S1, nR_S2] = trials2counts(data.stim.motionID - 1, data.response, data.conf, 2, 1);    
% [nR_S1, nR_S2] = trials2counts(data.stim.motionID(low_PE) - 1, data.response(low_PE), data.conf(low_PE), 2, 1);
  [nR_S1, nR_S2] = trials2counts(data.stim.motionID(low_conflict) - 1, data.response(low_conflict), data.conf(low_conflict), 2, 1);
    fit(sub) = fit_meta_d_MLE(nR_S1, nR_S2);
end

% mean_conf1 = mean(conf_1); mean_conf2 = mean(conf_2);
% mean_metada = 1.3786; % for low conflict condition 

% alpha = [0.1 0.3 0.5];
% [PE3, PE6, lowPE3, lowPE6] = fit_dprime();
alpha = [0.1];

for sub = 1:length(subjectIDs)
    for i = 1:length(alpha)
        prating = [conf_1(sub) conf_2(sub)];
        TI_KML2015_1A_md_Cx_search_loop(alpha(i), fit(sub).meta_da, prating, all_lowPE01(sub), subjectIDs{sub})
    end
end


%%
stim = [];

subjects = 20;

for s = 1:20
%   [PE6:lowC_lowCoh PE6:lowC_hiCoh PE3:hiC_lowCoh PE3:hiC_hiCoh]
% hi coh low conflict, PE:NE 6:1
 % hi coh hi conflict, PE:NE 3:1
% low coh low conflict, PE:NE 6:1
% low coh hi conflict, PE:NE 3:1
%     stim(1,:) = [PE01(2,2), PE01(2,2)/6]; 
%     stim(2,:) = [PE01(4,2), PE01(4,2)/3];
%     stim(3,:) = [PE01(1,2), PE01(1,2)/6];
%     stim(4,:) = [PE01(3,2), PE01(3,2)/3];
    stim(1,:) = [PE01(s,2), PE01(s,2)/6]; 
    stim(2,:) = [PE01(s,4), PE01(s,4)/3];
    stim(3,:) = [PE01(s,1), PE01(s,1)/6];
    stim(4,:) = [PE01(s,3), PE01(s,3)/3];
    all_stim{s} = stim; 
end

% all_stim = stim;

parfor s = 1:subjects
%     tic
%     waitbar(sub/10,h)
%     [boldLR{s}, perf{s}, correctLR{s}, saveWLR{s}] = simBOLD_iters2('LR', 0.2, 100);
    [boldLR_x{s}, perf_x{s}, correctLR_x{s}, saveWLR_x{s}] = simBOLD_iters2(all_stim{s}, 0.2, 100,'x');
    [boldLR_delta{s}, perf_delta{s}, correctLR_delta{s}, saveWLR_delta{s}] = simBOLD_iters2(all_stim{s}, 0.2, 100,'delta');
    
    % prepare for AUC
    [bin_levels_x{s}, activity_x{s}, response_x{s}, answer_LR_x{s}] = LR_sim_output(perf_x{s}, boldLR_x{s}, correctLR_x{s});
    [bin_levels_delta{s}, activity_delta{s}, response_delta{s}, answer_LR_delta{s}] = LR_sim_output(perf_delta{s}, boldLR_delta{s}, correctLR_delta{s});
%     time(s) = toc;
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

%% getting behavioral performance stats and plotting them
% 
% stim = [ 0.008571, 0.001429; ... % hi coh low conflict, PE:NE 6:1
%          0.007500, 0.002500; ...   % hi coh hi conflict, PE:NE 3:1
%          0.006857, 0.001143; ... % low coh low conflict, PE:NE 6:1
%          0.006000, 0.002000 ]; % low coh hi conflict, PE:NE 3:1
 
% variables are set up so out the 400, each increment of 100 is a
 % different stim, so need to separate them by stim:
stim = [1 100; 101 200; 201 300; 301 400];
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
xlim([0.8, 2.8]);
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
xlim([0.8, 2.8]);
xlabel('d''')
ylabel('mean conf d model')
title(['data for n = ' num2str(subjects)])

% legend(['low conflict (PE/NE = ' num2str(PE_NE_low) ')'], ['high conflict (PE/NE = ' num2str(PE_NE_high) ')'],'Location','NorthWest');
legend('low PE', 'high PE')

%% OLD, using low PE instead of low conflict 
% clear


% subjectIDs = {'101','105','110','111','203','204','205','206','207','208',...
% '209','213','215','217','218','220','221','222','223','224'};
% % subjectIDs = {'101'};
% 
% for sub = 1:length(subjectIDs)
%     data = TN_fMRI_getSubjectData(subjectIDs{sub}, 'LR');
%     conf_1(sub) = mean(data.conf==1);
%     conf_2(sub) = mean(data.conf==2);
% 
%  % meta d'
%     low_PE = data.stim.PE_NE_ratio == 3;
% %     [nR_S1, nR_S2] = trials2counts(data.stim.motionID - 1, data.response, data.conf, 2, 1);    
%     [nR_S1, nR_S2] = trials2counts(data.stim.motionID(low_PE) - 1, data.response(low_PE), data.conf(low_PE), 2, 1);
%     fit(sub) = fit_meta_d_MLE(nR_S1, nR_S2);
% end
% 
% 
% % alpha = [0.1 0.3 0.5];
% % [PE3, PE6, lowPE3, lowPE6] = fit_dprime();
% alpha = [0.1667; 0.333];
% 
% for sub = 1:length(subjectIDs)
%     for i = 1:length(alpha)
%         prating = [conf_1(sub) conf_2(sub)];
%         TI_KML2015_1A_md_Cx_search_loop(alpha(i), fit(sub).meta_da, prating, lowPE6(sub,2), lowPE3(sub,2), subjectIDs{sub})
%     end
% end