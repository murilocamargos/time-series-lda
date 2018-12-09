function [labels, theta, beta, zd] = tslda(TimeSeries, NumTopics, varargin)
% TSLDA  Utiliza��o do modelo Aloca��o Latente de Dirichlet (LDA, do ingl�s
% Latent Dirichlet Allocation) para agrupamento de s�ries temporais.
%
% Entradas
% ========
%   Required  TimeSeries      <array <numeric>>
%     A s�rie temporal a ser agrupada
%   Required  NumTopics       <integer>
%     O n�mero de grupos que se espera encontrar na s�rie
%   Optional  FuzzyForm       <string>
%     A forma fuzzy com que se deseja fuzzificar a s�rie temporal
%     Options: gauss (default), tri
%   Optional  NumWords        <integer>
%     O n�mero de palavras ou regras fuzzy que se deseja obter
%     Default: 5
%   Optional  NumWordsPerDoc  <integer>
%     O n�mero de palavras que se ter� em cada documento
%     Default: 10
%   Optional  Iterations      <integer>
%     O n�mero de itera��es para o amostrador de Gibbs
%     Default: 300
%   Optional  Alpha           <array <numeric>>
%     Um vetor com valores de alpha para cada t�pico ou um valor �nico que
%     ser� replicado para todos os t�picos
%     Default: 0.1
%   Optional  Gamma           <array <numeric>>
%     Um vetor com valores de gamma para cada palavra ou um valor �nico que
%     ser� replicado para todos as palavras
%     Default: 0.1
%
% Sa�das
% ======
%   labels <array <integer>>
%     R�tulos encontrados para cada valor da s�rie
%   theta  <array <numeric>>
%     Thetas encontrados pelo amostrador de Gibbs em cada instante
%   beta   <array <numeric>>
%     Betas encontrados pelo amostrador de Gibbs em cada instante
%   zd     <array <integer>>
%     Zds encontrados pelo amostrador de Gibbs em cada instante

    %% Valida��o dos par�metros do modelo
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

    %% Tranforma��o da s�rie temporal
    % Fuzzifica s�rie de entrada formando o universo de discurso
    words = fuzzify(p.Results.TimeSeries, p.Results.FuzzyForm, ...
        p.Results.NumWords);

    % Divide palavras do dicion�rio em documentos de tamanho igual e conta
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

    %% Inicializa��o dos par�metros do modelo LDA
    % O vetor de thetas � dado pela amostra de uma Dirichlet com parametros
    % alpha
    theta = zeros(p.Results.NumTopics, p.Results.Iterations + 1);
    theta(:, 1) = dirrnd(alpha, 1);

    % A matriz de betas � dado por amostras de uma Dirichlet com parametros
    % gamma
    beta = zeros(p.Results.NumWords, p.Results.NumTopics, ...
        p.Results.Iterations + 1);
    beta(:, 1, 1) = dirrnd(gamma, 1);
    for i = 2:p.Results.NumTopics
        beta(:, i, 1) = dirrnd(gamma, 1);
    end

    % O vetor inicial dos t�picos dos documentos � dado pela distribui��o
    % multinomial cujo par�metro � uma distribui��o uniforma, inicialmente.
    zd = zeros(numberOfDocs, p.Results.Iterations+1);
    zd(:, 1) = randsample(p.Results.NumTopics, numberOfDocs, true);

    %% Amostrados de Gibbs
    for i = 2:p.Results.Iterations+1
        % O novo vetor de alphas que servir� como par�metro para a nova
        % amostragem dos thetas depende da propor��o de t�picos atuais.
        newAlpha = alpha;
        topicCount = tabulate(zd(:, i-1));
        newAlpha(topicCount(:, 1)) = newAlpha(topicCount(:, 1)) + ...
            topicCount(:, 2)';
        theta(:, i) = dirrnd(newAlpha, 1);

        % Computa um novo gamma para cada t�pico com base na ocorr�ncia das
        % palavras em cada t�pico
        for k = 1:p.Results.NumTopics
            newGamma = sum(wordCount(zd(:, i-1) == k, :)) + gamma;
            beta(:, k, i) = dirrnd(newGamma, 1);
        end

        % Reamostra os t�picos Zd com base nos vetores de thetas e betas
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
%% Esta fun��o amostra de uma distribui��o de Dirichlet utilizando o
 % amostrador da distribui��o gama.
 
    k = length(thetas);
    p = zeros(n, k);
    for i = 1:k
        p(:, i) = gamrnd(thetas(i), 1, n, 1);
    end
    d = p./repmat(sum(p, 2), 1, k);
end

function w = fuzzify(y, form, num)
%% Esta fun��o fuzzifica a s�rie temporal criando um universo de discurso
 % em que as `num` regras fuzzy de forma `form` igualmente espa�adas s�o
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
