function [labels, theta, beta, zd] = tslda(TimeSeries, NumTopics, varargin)
%% TSLDA Uses the Latent Dirichlet Allocation (LDA) to grupo a fuzzy time
% series.
%
% Inputs
% ======
%   Required  TimeSeries      <array <numeric>>
%     The time-series to be clustered.
%   Required  NumTopics       <integer>
%     The number of topics you expect to find.
%   Optional  FuzzyForm       <string>
%     The fuzzy form used to fuzzify the time-series.
%     Options: gauss (default), tri
%   Optional  NumWords        <integer>
%     The numer of fuzzy words you wish to use.
%     Default: 5
%   Optional  NumWordsPerDoc  <integer>
%     The number of words per document.
%     Default: 10
%   Optional  Iterations      <integer>
%     Maximum number of iterations of the Gibbs sampler.
%     Default: 300
%   Optional  Alpha           <array <numeric>>
%     A vector with alpha values for each topic or a single value that
%     will be replicated for all topics.
%     Default: 0.1
%   Optional  Gamma           <array <numeric>>
%     A vector with gammas values for each topic or a single value that
%     will be replicated for all topics.
%     Default: 0.1
%
% Outputs
% =======
%   labels <array <integer>>
%     Labels found for each measurement.
%   theta  <array <numeric>>
%     Estimated thetas for each instant.
%   beta   <array <numeric>>
%     Estimated betas for each instant.
%   zd     <array <integer>>
%     Estimated zds for each instant.

    %% Parameter validation
    validPositiveInteger = @(x) isnumeric(x) && mod(x, 1) == 0 && (x > 0);
    validPositiveNumeric = @(x) isnumeric(x) && all(x > 0);
    validFuzzyForms = @(x) any(strcmp({'gauss', 'tri'}, x));

    p = inputParser;
    addRequired(p, 'TimeSeries', @isnumeric);
    addRequired(p, 'NumTopics', validPositiveInteger);
    addOptional(p, 'FuzzyForm', 'gauss', validFuzzyForms);
    addOptional(p, 'NumWords', 5, validPositiveInteger);
    addOptional(p, 'NumWordsPerDoc', 10, validPositiveInteger);
    addOptional(p, 'Iterations', 300, validPositiveInteger);
    addOptional(p, 'Alpha', 0.1, validPositiveNumeric);
    addOptional(p, 'Gamma', 0.1, validPositiveNumeric);
    parse(p, TimeSeries, NumTopics, varargin{:});

    % Validates the initial vector of alphas if provided
    alpha = p.Results.Alpha;
    if length(alpha) == 1
        alpha = ones(1, p.Results.NumTopics) * alpha;
    elseif size(p.Results.Alpha, 2) ~= p.Results.NumTopics
        error('The alpha vector should have the same size of the amount of topics.')
    end

    % Validates the initial vector of gammas if provided
    gamma = p.Results.Gamma;
    if length(gamma) == 1
        gamma = ones(1, p.Results.NumWords) * gamma;
    elseif size(p.Results.Gamma, 2) ~= p.Results.NumWords
        error('The gamma vector should have the same size of the amount of words.')
    end

    %% Time-series transformation
    % Fuzzyfies the input series
    words = fuzzify(p.Results.TimeSeries, p.Results.FuzzyForm, ...
        p.Results.NumWords);

    % Divides the dictionary words in document of equal size and counts
    % how mny words there is in each one of them
    numberOfDocs = ceil(length(words)/p.Results.NumWordsPerDoc);
    wordCount = zeros(numberOfDocs, p.Results.NumWords);
    for d = 1:numberOfDocs
        indexIni = (d - 1) * p.Results.NumWordsPerDoc + 1;
        indexEnd = min(indexIni + p.Results.NumWordsPerDoc - 1, ...
            length(p.Results.TimeSeries));
        wordsInDocD = tabulate(words(indexIni:indexEnd));
        wordCount(d, wordsInDocD(:, 1)) = wordsInDocD(:, 2);
    end

    %% Initializes the LDA parameters
    % The theta vector is given by a sample of the Dirichlet distribution
    % with parameters alpha
    theta = zeros(p.Results.NumTopics, p.Results.Iterations + 1);
    theta(:, 1) = dirrnd(alpha, 1);

    % The beta matrix is given by a sample of the Dirichlet distribution
    % with parameters gamma
    beta = zeros(p.Results.NumWords, p.Results.NumTopics, ...
        p.Results.Iterations + 1);
    beta(:, 1, 1) = dirrnd(gamma, 1);
    for i = 2:p.Results.NumTopics
        beta(:, i, 1) = dirrnd(gamma, 1);
    end

    % The initial vector of document topics is given by the multinomial
    % distribution with parameters following a uniform distribution.
    zd = zeros(numberOfDocs, p.Results.Iterations+1);
    zd(:, 1) = randsample(p.Results.NumTopics, numberOfDocs, true);

    %% Gibbs sampler
    for i = 2:p.Results.Iterations+1
        % The new alpha vector that will be used as parameters to sample
        % the thetas depende on the current topic proportional
        newAlpha = alpha;
        topicCount = tabulate(zd(:, i-1));
        newAlpha(topicCount(:, 1)) = newAlpha(topicCount(:, 1)) + ...
            topicCount(:, 2)';
        theta(:, i) = dirrnd(newAlpha, 1);

        % Compute a new gamma for each topic using the word cound on each
        % topic
        for k = 1:p.Results.NumTopics
            newGamma = sum(wordCount(zd(:, i-1) == k, :)) + gamma;
            beta(:, k, i) = dirrnd(newGamma, 1);
        end

        % Resample the topics Zd using the theta and beta vectors
        for d = 1:numberOfDocs
            probs = zeros(1, p.Results.NumTopics);
            for k = 1:p.Results.NumTopics
                probs(k) = log(theta(k, i)) + wordCount(d, :) * ...
                    log(beta(:, k, i));
            end

            probs = exp(probs);
            probs = probs/sum(probs);
            zd(d, i) = randsample(p.Results.NumTopics, 1, true, probs);
        end
    end

    labels = p.Results.TimeSeries * 0;
    for i = 1:length(zd(:, end))
        indexIni = (i - 1) * p.Results.NumWordsPerDoc + 1;
        indexEnd = min(indexIni + p.Results.NumWordsPerDoc - 1, ...
            length(p.Results.TimeSeries));
        labels(indexIni:indexEnd) = zd(i, end);
    end
end

function d = dirrnd(thetas, n)
%% DIRRND This function samples from a Dirichlet distribution using the
% gamma distribution sampler.
 
    k = length(thetas);
    p = zeros(n, k);
    for i = 1:k
        p(:, i) = gamrnd(thetas(i), 1, n, 1);
    end
    d = p./repmat(sum(p, 2), 1, k);
end

function w = fuzzify(y, form, num)
%% FUZZIFY This functino fuzzifies a time-series where a fixed number of
% fuzzy rules equaly spaced are used.
 
    limits = [floor(min(y)), ceil(max(y))];
    centers = linspace(min(limits), max(limits), num)';
    center_size = mean(diff(centers));

    if strcmp('gauss', form)
        sigma = center_size/3;
        rules = [ones(num, 1) * sigma, centers];
    elseif strcmp('tri', form)
        rules = [centers - center_size, centers, centers + center_size];
    end
    
    mdegree = zeros(length(y), num);
    for i = 1:num
        if strcmp('gauss', form)
            mdegree(:, i) = gaussmf(y, rules(i, :));
        elseif strcmp('tri', form)
            mdegree(:, i) = trimf(y, rules(i, :));
        end
    end

    [~, idx] = sort(mdegree, 2);
    w = idx(:, end);
end
