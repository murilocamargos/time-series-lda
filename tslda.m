function [labels, theta, beta, zd] = tslda(TimeSeries, NumTopics, varargin)
% TSLDA  Utilização do modelo Alocação Latente de Dirichlet (LDA, do inglês
% Latent Dirichlet Allocation) para agrupamento de séries temporais.
%
% Entradas
% ========
%   Required  TimeSeries      <array <numeric>>
%     A série temporal a ser agrupada
%   Required  NumTopics       <integer>
%     O número de grupos que se espera encontrar na série
%   Optional  FuzzyForm       <string>
%     A forma fuzzy com que se deseja fuzzificar a série temporal
%     Options: gauss (default), tri
%   Optional  NumWords        <integer>
%     O número de palavras ou regras fuzzy que se deseja obter
%     Default: 5
%   Optional  NumWordsPerDoc  <integer>
%     O número de palavras que se terá em cada documento
%     Default: 10
%   Optional  Iterations      <integer>
%     O número de iterações para o amostrador de Gibbs
%     Default: 300
%   Optional  Alpha           <array <numeric>>
%     Um vetor com valores de alpha para cada tópico ou um valor único que
%     será replicado para todos os tópicos
%     Default: 0.1
%   Optional  Gamma           <array <numeric>>
%     Um vetor com valores de gamma para cada palavra ou um valor único que
%     será replicado para todos as palavras
%     Default: 0.1
%
% Saídas
% ======
%   labels <array <integer>>
%     Rótulos encontrados para cada valor da série
%   theta  <array <numeric>>
%     Thetas encontrados pelo amostrador de Gibbs em cada instante
%   beta   <array <numeric>>
%     Betas encontrados pelo amostrador de Gibbs em cada instante
%   zd     <array <integer>>
%     Zds encontrados pelo amostrador de Gibbs em cada instante

    %% Validação dos parâmetros do modelo
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

    % Valida o vetor de inicial de alphas caso seja fornecido
    alpha = p.Results.Alpha;
    if length(alpha) == 1
        alpha = ones(1, p.Results.NumTopics) * alpha;
    elseif size(p.Results.Alpha, 2) ~= p.Results.NumTopics
        error('The alpha vector should have the same size of the amount of topics.')
    end

    % Valida o vetor de inicial de gammas caso seja fornecido
    gamma = p.Results.Gamma;
    if length(gamma) == 1
        gamma = ones(1, p.Results.NumWords) * gamma;
    elseif size(p.Results.Gamma, 2) ~= p.Results.NumWords
        error('The gamma vector should have the same size of the amount of words.')
    end

    %% Tranformação da série temporal
    % Fuzzifica série de entrada formando o universo de discurso
    words = fuzzify(p.Results.TimeSeries, p.Results.FuzzyForm, ...
        p.Results.NumWords);

    % Divide palavras do dicionário em documentos de tamanho igual e conta
    % quantas palavras de cada existem em cada documento
    numberOfDocs = ceil(length(words)/p.Results.NumWordsPerDoc);
    wordCount = zeros(numberOfDocs, p.Results.NumWords);
    for d = 1:numberOfDocs
        indexIni = (d - 1) * p.Results.NumWordsPerDoc + 1;
        indexEnd = min(indexIni + p.Results.NumWordsPerDoc - 1, ...
            length(p.Results.TimeSeries));
        wordsInDocD = tabulate(words(indexIni:indexEnd));
        wordCount(d, wordsInDocD(:, 1)) = wordsInDocD(:, 2);
    end

    %% Inicialização dos parâmetros do modelo LDA
    % O vetor de thetas é dado pela amostra de uma Dirichlet com parametros
    % alpha
    theta = zeros(p.Results.NumTopics, p.Results.Iterations + 1);
    theta(:, 1) = dirrnd(alpha, 1);

    % A matriz de betas é dado por amostras de uma Dirichlet com parametros
    % gamma
    beta = zeros(p.Results.NumWords, p.Results.NumTopics, ...
        p.Results.Iterations + 1);
    beta(:, 1, 1) = dirrnd(gamma, 1);
    for i = 2:p.Results.NumTopics
        beta(:, i, 1) = dirrnd(gamma, 1);
    end

    % O vetor inicial dos tópicos dos documentos é dado pela distribuição
    % multinomial cujo parâmetro é uma distribuição uniforma, inicialmente.
    zd = zeros(numberOfDocs, p.Results.Iterations+1);
    zd(:, 1) = randsample(p.Results.NumTopics, numberOfDocs, true);

    %% Amostrados de Gibbs
    for i = 2:p.Results.Iterations+1
        % O novo vetor de alphas que servirá como parâmetro para a nova
        % amostragem dos thetas depende da proporção de tópicos atuais.
        newAlpha = alpha;
        topicCount = tabulate(zd(:, i-1));
        newAlpha(topicCount(:, 1)) = newAlpha(topicCount(:, 1)) + ...
            topicCount(:, 2)';
        theta(:, i) = dirrnd(newAlpha, 1);

        % Computa um novo gamma para cada tópico com base na ocorrência das
        % palavras em cada tópico
        for k = 1:p.Results.NumTopics
            newGamma = sum(wordCount(zd(:, i-1) == k, :)) + gamma;
            beta(:, k, i) = dirrnd(newGamma, 1);
        end

        % Reamostra os tópicos Zd com base nos vetores de thetas e betas
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
%% Esta função amostra de uma distribuição de Dirichlet utilizando o
 % amostrador da distribuição gama.
 
    k = length(thetas);
    p = zeros(n, k);
    for i = 1:k
        p(:, i) = gamrnd(thetas(i), 1, n, 1);
    end
    d = p./repmat(sum(p, 2), 1, k);
end

function w = fuzzify(y, form, num)
%% Esta função fuzzifica a série temporal criando um universo de discurso
 % em que as `num` regras fuzzy de forma `form` igualmente espaçadas são
 % utilizadas.
 
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
