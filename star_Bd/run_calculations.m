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
fpC = zeros(1, length(selection_intensities));
% fixation probabilities of non-producers
fpD = zeros(1, length(selection_intensities));

% loop through selection intensities. change 'parfor' to 'for' for serial
% loop.
for i = 1:length(selection_intensities)
    delta = selection_intensities(i);
    
    % initialize (empty) transition matrix
    transition_matrix = [];
    if strcmp(social_good, 'pp')
        transition_matrix = build_matrix_pp(N, b, c, delta);
    elseif strcmp(social_good, 'cf')
        transition_matrix = build_matrix_cf(N, b, c, delta);
    elseif strcmp(social_good, 'ff')
        transition_matrix = build_matrix_ff(N, b, c, delta);
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
        muC_H = 1/N;
        muC_L = 1-1/N;
        muD_H = 1/N;
        muD_L = 1-1/N;
    elseif strcmp(mutant_appearance, 'temperature')
        muC_H = 1-1/N;
        muC_L = 1/N;
        if strcmp(social_good, 'pp')
            pay_H = (N-1)*(b-c);
            pay_L = b-c;
            fec_H = exp(delta*(pay_H-max(pay_H, pay_L)));
            fec_L = exp(delta*(pay_L-max(pay_H, pay_L)));
            muD_H = (N-1)*fec_L/(fec_H+(N-1)*fec_L);
            muD_L = fec_H/(fec_H+(N-1)*fec_L);
        elseif strcmp(social_good, 'cf')
            pay_H = (N-1)*b-c;
            pay_L = -c;
            pay_L_recipient = b-c;
            fec_H = exp(delta*(pay_H-max([pay_H, pay_L, pay_L_recipient])));
            fec_L = exp(delta*(pay_L-max([pay_H, pay_L, pay_L_recipient])));
            fec_L_recipient = exp(delta*(pay_L_recipient-max([pay_H, pay_L, pay_L_recipient])));
            muD_H = ((N-2)*fec_L+fec_L_recipient)/(fec_H+(N-2)*fec_L+fec_L_recipient);
            muD_L = fec_H/(fec_H+(N-2)*fec_L+fec_L_recipient);
        elseif strcmp(social_good, 'ff')
            pay_H = (N-1)*b-c;
            pay_L = b/(N-1)-c;
            fec_H = exp(delta*(pay_H-max(pay_H, pay_L)));
            fec_L = exp(delta*(pay_L-max(pay_H, pay_L)));
            muD_H = (N-1)*fec_L/(fec_H+(N-1)*fec_L);
            muD_L = fec_H/(fec_H+(N-1)*fec_L);
        else
            error('Unrecognized social good.');
        end
    else
        error('Unrecognized appearance distribution.');
    end
    
    % calculate expected fixation probabilities
    fpC(i) = muC_L*fpC_leaf + muC_H*fpC_hub;
    fpD(i) = muD_L*fpD_leaf + muD_H*fpD_hub;
end

% plot results
plot(selection_intensities, fpD, 'Color', [0.6, 0, 0], 'LineWidth', 2);
hold on;
plot(selection_intensities, fpC, 'Color', [0, 0, 0.6], 'LineWidth', 2);
axis square; box on; grid on;
