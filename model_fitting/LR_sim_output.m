function [bin_levels, activity, response, answer_LR] = LR_sim_output(perf, boldLR, correctLR)

% bin_levels, activity] = LR_sim_output(perf, boldLR);
% takes in behavioral performance and simulated bold activity and outputs
% binarized confidence levels and activity to be of the right shape for AUC
% classification

confidenceX = [];
confidenceD = [];
response = [];
activity = [];
answer_LR = [];

for i = 1:(size(perf,1))
    for j = 1:(size(perf,2))
        response = [response;perf(i,j).resp];
        confidenceX = [confidenceX;perf(i,j).conf_Cx];
        confidenceD = [confidenceD;perf(i,j).conf_Cd];
        activity = [activity;(boldLR(i,:,j))]; %becomes stacked in a way where 
        % 1:1000 is stim comb1, 1001:2000 is stim comb2, etc.
        answer_LR = [answer_LR;correctLR(i,j)];
    end
end

% x model
q = quantile(confidenceX(:),3); 
levelsX = zeros(size(confidenceX));
levelsX(confidenceX < q(1) & confidenceX >= 0) = 1;
levelsX(confidenceX < q(2) & confidenceX >= q(1)) = 2;
levelsX(confidenceX < q(3) & confidenceX >= q(2)) = 3;
levelsX(confidenceX >= q(3)) = 4;
bin_levelsX = levelsX > 2;

% delta model
q = quantile(confidenceD(:),3); 
levelsD = zeros(size(confidenceD));
levelsD(confidenceD < q(1) & confidenceD >= 0) = 1;
levelsD(confidenceD < q(2) & confidenceD >= q(1)) = 2;
levelsD(confidenceD < q(3) & confidenceD >= q(2)) = 3;
levelsD(confidenceD >= q(3)) = 4;
bin_levelsD = levelsD > 2;

bin_levels.X = bin_levelsX;
bin_levels.D = bin_levelsD;

%show response correct as a function of condition
% mean_resp = mean(response==0,2);
% mean_confX = mean(levelsX,2);
% mean_confD = mean(levelsD,2);
% mean_activity = mean(boldLR,3);
end