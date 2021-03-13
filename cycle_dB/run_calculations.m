% RUN_CALCULATIONS
%   main script to run calculations on the cycle for dB updating.
%   parameters can be changed as needed.

addpath('../');

% population size
N = 10;

% state space for the process
states = transpose(0:1:N);

% benefit and cost for producers
b = 4;
c = 1;

% range of selection intensities
selection_intensities = 0.001:0.001:3.999;

% type of social good (either 'pp', 'cf', or 'ff')
social_good = 'ff';

% fixation probabilities of producers
fixation_probabilities_C = zeros(1, length(selection_intensities));
% fixation probabilities of non-producers
fixation_probabilities_D = zeros(1, length(selection_intensities));

% loop through selection intensities. change 'parfor' to 'for' for serial
% loop.
parfor i = 1:length(selection_intensities)
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
    fixation_probabilities_C(i) = fixation_probability(1, N, states, transition_matrix);
    fixation_probabilities_D(i) = fixation_probability(N-1, 0, states, transition_matrix);
end

% plot results
plot(selection_intensities, fixation_probabilities_D, 'Color', [0.6, 0, 0], 'LineWidth', 2);
hold on;
plot(selection_intensities, fixation_probabilities_C, 'Color', [0, 0, 0.6], 'LineWidth', 2);
axis square; box on; grid on;
