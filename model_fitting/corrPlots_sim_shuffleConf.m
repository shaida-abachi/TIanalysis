function [ttest_results, paired_ttest, corrs] = corrPlots_sim_shuffleConf(subjects, dolog, makePlots)

% corrPlots plots the correlation and slope between inhibitory index and AUC
% subjects: a string of subjects of length 1 x n: e.g. {'100','101'}
% dolog: either 1 or 0 to take the log of inhibitory index

% cd auc
%close all; clc
if nargin < 1
     cd auc/
    subjects = {'101','105','110','111','203','204','205','206','207',...
        '208','209','213','215','217','218','220','221','222','223','224'}; 
    dolog=1; makePlots = 0;
end

%below for test 
% clear all; close all; clc
% subjects = {'101','105','110','111','203','204','205','206','207','208',...
%     '209','213','215','217','218','220','221','222','223','224'}; ...
%      dolog=1; makePlots = 0;

for s = 1:length(subjects)
    sub = subjects{s};
%     fnames{s} = ['Iinx_AUC_' sub '_minusoverplus_20densitythresh.mat'];
    fnames{s} = ['Iinx_AUC_' sub '_minusoverplus.mat']; % 25 density thresh
end

symbols = {'x','o','k','s','d','p','h','*','+','>','<','v','_','|','^','x','o','k','s','d',};
symbols = symbols(1,1:length(subjects));

%%

for s = 1:length(subjects)
%     load data
    load(fnames{s});
    AUCconf = AUCconf(randperm(length(AUCconf)));
    AUCdec;
    IDX = inh_inx;
    %IDX;

        % do we take the log?
    if dolog, inhIndx = log10(IDX); else, inhIndx = IDX; end

      all_conf_AUCs{s} = AUCconf';
      conf_means(s) = mean(AUCconf');
      
      all_dec_AUCs{s} = AUCdec';
      dec_means(s) = mean(AUCdec');
      all_inh_inx{s} = inhIndx';
    
     conf_minus_deci = AUCconf - AUCdec; %NEW 04202023
     all_conf_minus_deci{s} = conf_minus_deci';
      f_minus = fit(inhIndx',conf_minus_deci','poly1');
      fit_minus_vals.c(s,:) = coeffvalues(f_minus);
      fits_minus.c{s,:} = f_minus;
      
      if makePlots
          figure(40)
          h = plot(f_minus,inhIndx',conf_minus_deci');
          set(h,'LineWidth',1)
          title('conf - deci')
          xlabel('log inh index')
          ylabel('AUC')
          hold on
      end
    
    
      %dense
      [r,p] = corr(inhIndx',conf_minus_deci',"type",'Spearman');
%       disp(['S=' num2str(s) ', R_conf=' num2str(r) ', p=' num2str(p)])
      spearMinus.c(s,:) = [r,p];
    
      [r,p] = corr(inhIndx',conf_minus_deci',"type",'Pearson');
%       disp(['S=' num2str(s) ', R_conf=' num2str(r) ', p=' num2str(p)])
      pearsMinus.c(s,:) = [r,p];
    
    [~,conf_p(s)] = ttest(AUCconf,0.5);
    conf_means(s) = mean(AUCconf);
    [~,dec_p(s)] = ttest(AUCdec,0.5);
    dec_means(s) = mean(AUCdec);

    % plot InhIndex vs AUC for decision vs conf
      if makePlots
          figure(1)
          subplot(1,2,1)
          hold on
          scatter(inhIndx,AUCconf,symbols{s})
          title('Confidence')
          hold on
          yline(0.5,'r--', 'LineWidth', 1);
      end

      f = fit(inhIndx',AUCconf','poly1');
      fit_vals.c(s,:) = coeffvalues(f);
      fits.c{s,:} = f;

      if makePlots
          h = plot(f,inhIndx',AUCconf');
          set(h,'LineWidth',1)
          if dolog, xlabel('log InhIndex'); else, xlabel('InhIndex'); end
          ylabel('ConfAUC')
          
          if s == length(subjects)
              CONFIDENCE = cat( 1, all_conf_AUCs{:} );
              INHIB_INX = cat(1,all_inh_inx{:});
              c = fit(INHIB_INX,CONFIDENCE,'poly1');
              cp = plot(c,INHIB_INX,CONFIDENCE);
              set(cp,'LineWidth',4)
          end
      end

      %hold off
      
      [r,p] = corr(inhIndx',AUCconf',"type",'Spearman');
%       disp(['S=' num2str(s) ', R_conf=' num2str(r) ', p=' num2str(p)])
      spear.c(s,:) = [r,p];

      [r,p] = corr(inhIndx',AUCconf',"type",'Pearson');
%       disp(['S=' num2str(s) ', R_conf=' num2str(r) ', p=' num2str(p)])
      pears.c(s,:) = [r,p];
%       makePretty(20)

      if makePlots
          subplot(1,2,2)
          hold on
          scatter(inhIndx,AUCdec,symbols{s})
          title('Decision')
          hold on
          yline(0.5,'r--', 'LineWidth', 1);
      end

      fd = fit(inhIndx',AUCdec','poly1');
      fit_vals.d(s,:) = coeffvalues(fd);
      fits.d{s,:} = fd;
      
      if makePlots
          h = plot(fd,inhIndx',AUCdec');
          set(h,'LineWidth',1)
          if dolog, xlabel('log InhIndex'); else, xlabel('InhIndex'); end
          ylabel('DecAUC') 
      end

      [r,p] = corr(inhIndx',AUCdec',"type",'Spearman');
%       disp(['S=' num2str(s) ', R_dec=' num2str(r) ', p=' num2str(p)])
      spear.d(s,:) = [r,p];

      [r,p] = corr(inhIndx',AUCdec',"type",'Pearson');
%       disp(['S=' num2str(s) ', R_dec=' num2str(r) ', p=' num2str(p)])
      pears.d(s,:) = [r,p];

      if makePlots
          if s == length(subjects)
              DECISION = cat( 1, all_dec_AUCs{:} );
              d = fit(INHIB_INX,DECISION,'poly1');
              dp = plot(d,INHIB_INX,CONFIDENCE);
              set(dp,'LineWidth',4)
    
              figure(41)
              ALL = cat(1,all_conf_minus_deci{:});
              INHIB_INX = cat(1,all_inh_inx{:});
              c = fit(INHIB_INX,ALL,'poly1');
              cp = plot(c,INHIB_INX,ALL);
              set(cp,'LineWidth',4)
          end

          hold off
          makePretty(20)
      
      end

      [~,conf_deci_ttest(s)] = ttest(AUCconf',AUCdec');

end

corrs.spear = spear;
corrs.pears = pears;
corrs.fit_vals = fit_vals;

%%
spearTrans.c(:,1) = atanh(spear.c(:,1));
spearTrans.d(:,1) = atanh(spear.d(:,1));
pearsTrans.c(:,1) = atanh(pears.c(:,1));
pearsTrans.d(:,1) = atanh(pears.d(:,1));

if makePlots
    figure(2)
    hold on
    % bar([corrs.c(:,1) corrs.d(:,1)])
    bar([mean(spearTrans.c(:,1)),mean(spearTrans.d(:,1))])
    errorbar([mean(spearTrans.c(:,1)),mean(spearTrans.d(:,1))],[std(spearTrans.c(:,1)),mean(spearTrans.d(:,1))]./sqrt(length(subjects)),'k','LineStyle','none')
    set(gca,'xtick',1:2,'xticklabel',{'confidence','decision'})
    ylabel('Spearmans Rho')
    title('Spearman Corr')
    makePretty(20)
end
[~,spearTransp2] = ttest(spearTrans.c(:,1));
[~,spearTransp3] = ttest(spearTrans.d(:,1));
[spearTrans_h,spearTrans_p] = ttest(spearTrans.c(:,1),spearTrans.d(:,1));

if makePlots
    figure(3)
    hold on
    % bar([corrs.c(:,1) corrs.d(:,1)])
    bar([mean(pearsTrans.c(:,1)),mean(pearsTrans.d(:,1))])
    errorbar([mean(pearsTrans.c(:,1)),mean(pearsTrans.d(:,1))],[std(pearsTrans.c(:,1)),mean(pearsTrans.d(:,1))]./sqrt(length(subjects)),'k','LineStyle','none')
    set(gca,'xtick',1:2,'xticklabel',{'confidence','decision'})
    ylabel('Pearsons R')
    title('Pearson Corr')
    makePretty(20)
end
% 
% z.c = atanh(pears.c(:,1));
% z.d = atanh(pears.d(:,1));
[~,pearsTransp2] = ttest(pearsTrans.c(:,1));
[~,pearsTransp3] = ttest(pearsTrans.d(:,1));
%paired samples ttest
[pearsTrans_h,pearsTrans_p] = ttest(pearsTrans.c(:,1),pearsTrans.d(:,1));

if makePlots
    figure(4)
    hold on
    bar([mean(fit_vals.c(:,1)),mean(fit_vals.d(:,1))])
    errorbar([mean(fit_vals.c(:,1)),mean(fit_vals.d(:,1))],[std(fit_vals.c(:,1)),mean(fit_vals.d(:,1))]./sqrt(length(subjects)),'k','LineStyle','none')
    set(gca,'xtick',1:2,'xticklabel',{'confidence','decision'})
    ylabel('fitted slope')
    title('Linear fit slope')
    makePretty(20)
end

[~,fitp2] = ttest(fit_vals.c(:,1));
[~,fitp3] = ttest(fit_vals.d(:,1));

ttest_results = [spearTransp2;spearTransp3;pearsTransp2;pearsTransp3;fitp2;fitp3];

%paired samples ttest
[fit_h,fit_p] = ttest(fit_vals.c(:,1),fit_vals.d(:,1));

paired_ttest = [spearTrans_p;pearsTrans_p;fit_p];

% display('hi')

end