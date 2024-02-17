function [bold, savePerf, correctAns, saveW] = simBOLD_iters2(stim, inh, iters, model, varargin)

% simBOLD_iters2 runs through simBOLD for many iterations
% [bold, perf, saveW] = simBOLD_iters(stim, inh, iters, behavior, varargin)
% INPUTS-
%   stim: can either be 'MI', 'LR', or a ~ x 2 vector of PE/NE fractions
%   inh: inhibitory interneuron percentages. can be a single number or a
%       vector of many inhib interneuron percentages
%   iters: number of iterations
%   behavior: whether behavioral response is output
%
% OUTPUTS-
%   bold: simulated BOLD activity that is size of columns of stim x inh x
%   (v vs delta percentage) x iters
%   perf: behavioral performance and confidence of simulated subject

if ~exist('stim')
    stim = [0.01, 0];
end

if ~exist('inh') 
    inh = 0.2;
end

if ~exist('iters')
    iters = 1;
end

if ~exist('model')
    model = 'x';
end

if stim == 'MI'
    stim = [ 0.01, 0;      ... % high energy, low conflict
                0.007, 0.003; ... % high energy, high conflict
                0.007, 0;     ... % low energy, low conflict
                0.0042, 0.0028];  % low energy, high conflict
elseif stim == 'LR' %
    stim = [ 0.008571, 0.001429; ... % hi coh low conflict, PE:NE 6:1
                0.007500, 0.002500; ...   % hi coh hi conflict, PE:NE 3:1
                0.006857, 0.001143; ... % low coh low conflict, PE:NE 6:1
                0.006000, 0.002000]; % low coh hi conflict, PE:NE 3:1
else 
    stim = stim;
end

% if model == 'x' SPLA commented out 072023 because TI_sim_trial now does
% both
%     param.tau = 52;
% elseif model == 'delta'
%     param.tau = 21;
% end

stim_PE_NEs = stim; 
%stim_PE_NEs = stim_PE_NEs * 10;
           
x_vs_delta_freqs = linspace(1, 0, 20);  % relative frequency of x vs delta neurons in a voxel, when considering only non-inhibitory neurons
inh_int_fractions = inh;
n_inh    = length(inh_int_fractions);
n_stim   = size(stim_PE_NEs,1);
n_x_vs_d = length(x_vs_delta_freqs);

% loop over the factors

saveW = [];

param.tau   = 52; %20; SPLA edited 072023
param.T     = 1;   % decision threshold
param.sigma = 0.1; % standard deviation of accumulation noise
param.tmax  = 1e4; % max # of timesteps before simulation exits
%param.S_i = stim; %above taken from TI_KML2015_1A_d_search_loop, S_i taken from example
param.save_traces = 1;
param.exit_on_response = 0;
%[perf, traces] = TI_sim_trial(param); % DONT CALL THIS HERE MAYBE?

%h = waitbar(0,'Simulating BOLD...');
tic
for i_trial = 1:iters
%     display(i_trial)
    %waitbar(i_trial/iters,h)
    for i_stim = 1:n_stim
%         display(i_stim)
        if rand > .5
            param.S_i = stim_PE_NEs(i_stim,:)';
            correctAns(i_stim, i_trial) = 0;
        else
            param.S_i = flipud(stim_PE_NEs(i_stim,:)');
            correctAns(i_stim, i_trial) = 1;
        end
        [perf, traces] = TI_sim_trial(param); 
        savePerf(i_stim, i_trial) = perf;

        for i_x = 1:n_x_vs_d
%             display(i_x)
            % get stim/inh/x_vs_d settings for current simulation
            inh_int_fraction = 0.2;
            x_vs_delta_freq  = x_vs_delta_freqs(i_x);
            
            % figure out relative fraction of x and delta neurons
            x_fraction     = (1 - inh_int_fraction) * x_vs_delta_freq;
            delta_fraction = 1 - inh_int_fraction - x_fraction;
            
            %w = [0, x_fraction, inh_int_fraction, delta_fraction];
            W.Ws = 0;
            W.Wx = x_fraction;
            W.Wy = inh_int_fraction;
            W.Wdelta = delta_fraction;
            saveW = [saveW;W];

            % sim bold
            bold(i_stim, i_x, i_trial) = simBOLD_only(traces, W);
        end
    end
end
toc
end
