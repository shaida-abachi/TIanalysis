function [] = plotLRbehavior(subjects, response, answer_LR, bin_levels, changeLims, figTxtSize, varargin)
% Plot behavioral performance stats for low and high PE conditions
% subjects: number of simulated subjects
% response, answer_LR, bin_lebels all output from simBOLDiters
% changeLims: true or false, to specify limits of plots
% figTxtSize: size of text in figures, using makePretty.m function


% stim = [ 0.008571, 0.001429; ... % hi coh low conflict, PE:NE 6:1
%          0.007500, 0.002500; ...   % hi coh hi conflict, PE:NE 3:1
%          0.006857, 0.001143; ... % low coh low conflict, PE:NE 6:1
%          0.006000, 0.002000 ]; % low coh hi conflict, PE:NE 3:1
 
if ~exist('changeLims')
    changeLims = 1; % specify limits of plots?
end

if ~exist('figTxtSize')
        figTxtSize = 15; % size of text in figures
end

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

subplot(1,3,1);
% bar_error(d_m2, d_sem2);
barWithErrors(d_m2, d_sem2);
if changeLims
    set(gca, 'ylim', [0.8 2.8]);
    set(gca,'XTick',[1 2])
end
% set(gca,'XTickLabel',{'low conflict','hi conflict'})
set(gca,'XTickLabel',{'high PE:NE','low PE:NE'})
ylabel('d''')
% title(['data for n = ' num2str(subjects)])

subplot(1,3,2);
% bar_error(mc_m_x2, mc_sem_x2);
barWithErrors(mc_m_x2, mc_sem_x2);
yy = get(gca, 'ylim');
% set(gca, 'ylim', [1 yy(2)]);
if changeLims
    set(gca, 'ylim', [1.3 1.7]);
    set(gca,'XTick',[1 2])
end
% set(gca,'XTickLabel',{'low conflict','hi conflict'})
set(gca,'XTickLabel',{'high PE:NE','low PE:NE'})
legend('low coherence','high coherence')
ylabel('mean conf x')

subplot(1,3,3);
% bar_error(mc_m_d2, mc_sem_d2);
barWithErrors(mc_m_d2, mc_sem_d2);
yy = get(gca, 'ylim');
% set(gca, 'ylim', [1 yy(2)]);
if changeLims
    set(gca, 'ylim', [1.3 1.7]);
    set(gca,'XTick',[1 2])
end
% set(gca,'XTickLabel',{'low conflict','hi conflict'})
set(gca,'XTickLabel',{'high PE:NE','low PE:NE'})
legend('low coherence','high coherence')
ylabel('mean conf δ')
makePretty(figTxtSize)

% plot conf_x vs d'
figure; subplot(1,2,1); hold on;
plot(d_m2(1,:), mc_m_x2(1,:), 'go-', 'LineWidth', 2);
plot(d_m2(2,:), mc_m_x2(2,:), 'ro-', 'LineWidth', 2);
hv_errorbars(d_m2(1,:), mc_m_x2(1,:), d_sem2(1,:), mc_sem_x2(1,:), 'g-', .5);
hv_errorbars(d_m2(2,:), mc_m_x2(2,:), d_sem2(2,:), mc_sem_x2(2,:), 'r-', .5);
if changeLims
    xlim([0.8, 2.8]);
    ylim([1.3, 1.7]);
end
xlabel('d''')
ylabel('mean conf x')
% title(['data for n = ' num2str(subjects)])

% legend(['low conflict (PE/NE = ' num2str(PE_NE_low) ')'], ['high conflict (PE/NE = ' num2str(PE_NE_high) ')'],'Location','NorthWest');
% legend('low conflict', 'hi conflict')
legend('high PE:NE', 'low PE:NE')

% plot conf_d vs d'
% figure; hold on;
subplot(1,2,2); hold on;
plot(d_m2(1,:), mc_m_d2(1,:), 'go-', 'LineWidth', 2);
plot(d_m2(2,:), mc_m_d2(2,:), 'ro-', 'LineWidth', 2);
hv_errorbars(d_m2(1,:), mc_m_d2(1,:), d_sem2(1,:), mc_sem_d2(1,:), 'g-', .5);
hv_errorbars(d_m2(2,:), mc_m_d2(2,:), d_sem2(2,:), mc_sem_d2(2,:), 'r-', .5);
if changeLims
    xlim([0.8, 2.8]);
    ylim([1.3, 1.7]);
end
xlabel('d''')
ylabel('mean conf δ')
legend('high PE:NE', 'low PE:NE')
% title(['data for n = ' num2str(subjects)])
% legend(['low conflict (PE/NE = ' num2str(PE_NE_low) ')'], ['high conflict (PE/NE = ' num2str(PE_NE_high) ')'],'Location','NorthWest');
% legend('low conflict', 'hi conflict')

makePretty(figTxtSize)

end