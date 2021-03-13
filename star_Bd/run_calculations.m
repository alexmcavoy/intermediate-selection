% RUN_CALCULATIONS
%   main script to run calculations on the star for Bd updating.
%   parameters can be changed as needed.

addpath('../');

% population size
N = 10;

% state space for the process
states = cartesian_product(transpose(0:1:1), transpose(0:1:N-1));

% benefit and cost for producers
b = 5;
c = 1;

% range of selection intensities
selection_intensities = 0.001:0.001:1.999;

% type of social good (either 'pp', 'cf', or 'ff')
social_good = 'ff';

% type of mutant-appearance distribution (either 'uniform' or 'temperature')
mutant_appearance = 'uniform';

% fixation probabilities of producers
fixation_probabilities_C = zeros(1, length(selection_intensities));
% fixation probabilities of non-producers
fixation_probabilities_D = zeros(1, length(selection_intensities));

% loop through selection intensities. change 'parfor' to 'for' for serial
% loop.
for i = 1:length(selection_intensities)
    intensity = selection_intensities(i);
    
    % initialize (empty) transition matrix
    transition_matrix = [];
    if strcmp(social_good, 'pp')
        transition_matrix = build_matrix_pp(N, b, c, intensity);
    elseif strcmp(social_good, 'cf')
        transition_matrix = build_matrix_cf(N, b, c, intensity);
    elseif strcmp(social_good, 'ff')
        transition_matrix = build_matrix_ff(N, b, c, intensity);
    else
        error('Unrecognized social good.');
    end
    
    % probability that a single producer at the hub fixes
    initial_state = [1, 0];
    final_state = [1, N-1];
    fpC_hub = fixation_probability(initial_state, final_state, states, transition_matrix);

    % probability that a single non-producer at the hub fixes
    initial_state = [0, N-1];
    final_state = [0, 0];
    fpD_hub = fixation_probability(initial_state, final_state, states, transition_matrix);
    
    % probability that a single producer at a leaf fixes
    initial_state = [0, 1];
    final_state = [1, N-1];
    fpC_leaf = fixation_probability(initial_state, final_state, states, transition_matrix);

    % probability that a single non-producer at a leaf fixes
    initial_state = [1, N-2];
    final_state = [0, 0];
    fpD_leaf = fixation_probability(initial_state, final_state, states, transition_matrix);

    % calculate weights for the hub and leaf in the two monomorphic states
    if strcmp(mutant_appearance, 'uniform')
        muC_hub = 1/N;
        muC_leaf = 1-1/N;
        muD_hub = 1/N;
        muD_leaf = 1-1/N;
    elseif strcmp(mutant_appearance, 'temperature')
        muC_hub = 1-1/N;
        muC_leaf = 1/N;
        if strcmp(social_good, 'pp')
            payoff_hub = (N-1)*(b-c);
            payoff_leaf = b-c;
            fecundity_hub = exp(intensity*(payoff_hub-max(payoff_hub, payoff_leaf)));
            fecundity_leaf = exp(intensity*(payoff_leaf-max(payoff_hub, payoff_leaf)));
            muD_hub = (N-1)*fecundity_leaf/(fecundity_hub+(N-1)*fecundity_leaf);
            muD_leaf = fecundity_hub/(fecundity_hub+(N-1)*fecundity_leaf);
        elseif strcmp(social_good, 'cf')
            payoff_hub = (N-1)*b-c;
            payoff_leaf = -c;
            payoff_leaf_recipient = b-c;
            fecundity_hub = exp(intensity*(payoff_hub-max([payoff_hub, payoff_leaf, payoff_leaf_recipient])));
            fecundity_leaf = exp(intensity*(payoff_leaf-max([payoff_hub, payoff_leaf, payoff_leaf_recipient])));
            fecundity_leaf_recipient = exp(intensity*(payoff_leaf_recipient-max([payoff_hub, payoff_leaf, payoff_leaf_recipient])));
            muD_hub = ((N-2)*fecundity_leaf+fecundity_leaf_recipient)/(fecundity_hub+(N-2)*fecundity_leaf+fecundity_leaf_recipient);
            muD_leaf = fecundity_hub/(fecundity_hub+(N-2)*fecundity_leaf+fecundity_leaf_recipient);
        elseif strcmp(social_good, 'ff')
            payoff_hub = (N-1)*b-c;
            payoff_leaf = b/(N-1)-c;
            fecundity_hub = exp(intensity*(payoff_hub-max(payoff_hub, payoff_leaf)));
            fecundity_leaf = exp(intensity*(payoff_leaf-max(payoff_hub, payoff_leaf)));
            muD_hub = (N-1)*fecundity_leaf/(fecundity_hub+(N-1)*fecundity_leaf);
            muD_leaf = fecundity_hub/(fecundity_hub+(N-1)*fecundity_leaf);
        else
            error('Unrecognized social good.');
        end
    else
        error('Unrecognized appearance distribution.');
    end
    
    % calculate expected fixation probabilities
    fixation_probabilities_C(i) = muC_leaf*fpC_leaf + muC_hub*fpC_hub;
    fixation_probabilities_D(i) = muD_leaf*fpD_leaf + muD_hub*fpD_hub;
end

% plot results
plot(selection_intensities, fixation_probabilities_D, 'Color', [0.6, 0, 0], 'LineWidth', 2);
hold on;
plot(selection_intensities, fixation_probabilities_C, 'Color', [0, 0, 0.6], 'LineWidth', 2);
axis square; box on; grid on;
