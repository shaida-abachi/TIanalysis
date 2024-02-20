function [lowPE, highPE] = fit_dprime(low_high)
% fits stim (S) to mean d' of all conditions 
% low_high is a 1 x 2 vector that is either [3 6] or [5 10] corresponding
% to alphas of [0.333 0.1667] or [0.2 0.1]
% outputs fitted S of the 4 conditons for the low and high conditions 

% 0.1 alpha of 0.1 for high PE condition (1/10)
fit_1 = [-879.941967991654 136.383362056827	-0.0103557840141724];

% 0.2 alpha of 0.2 for low PE condition (1/5)
fit_2 = [-665.588899044654 119.317503767596 -0.00453928763253160];

% alpha of 0.1667 (1/6) for high PE condition
fit_1_6 = [-728.325283097094 125.445649494197 -0.00316806314821995];

% alpha of 0.333 (1/3) for low PE condition
fit_1_3 = [-447.471395516241 97.6413608630508 0.00332510875338433];

% d' of all 20 subjects in response to the 4 stimulus conditions 

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

 all_dp = mean(sub,1); % mean across subjects (output is 1 x 4)

for stim = 1:4
    p1 = fit_1;
    p1(3) = fit_1(3) -  all_dp(stim);
    PE01(stim,:) = roots(p1); 

    p2 = fit_2;
    p2(3) = fit_2(3) -  all_dp(stim);
    PE02(stim,:) = roots(p2);

    p1_6 = fit_1_6;
    p1_6(3) = fit_1_6(3) -  all_dp(stim);
    PE01_6(stim,:) = roots(p1_6);

    p1_3 = fit_1_3;
    p1_3(3) = fit_1_3(3) -  all_dp(stim);
    PE01_3(stim,:) = roots(p1_3); 
end

if low_high == [3 6]
    lowPE = PE01_3;
    highPE = PE01_6;
elseif low_high == [5 10]
    lowPE = PE02;
    highPE = PE01;
end

end
