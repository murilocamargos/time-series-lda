dirrnd <- function(thetas, n) {
  # Esta função amostra de uma distribuição de Dirichlet utilizando o
  # amostrador da distribuição gama.
  k <- length(thetas)
  p <- matrix(numeric(n*k), nrow=n)
  for (i in 1:k) {
    p[, i] = rgamma(n, thetas[i])
  }
  return(p/matrix(rep(rowSums(p), k), nrow=n))
}

fuzzify <- function(y, form, num) {
  # Esta função fuzzifica a série temporal criando um universo de discurso
  # em que as `num` regras fuzzy de forma `form` igualmente espaçadas são
  # utilizadas.
  
  limits <- c(floor(min(y)), ceiling(max(y)));
  centers <- seq(min(limits), max(limits), length=num)
  center_size <- mean(diff(centers))
  
  if (form == 'gauss') {
    sigma <- center_size/3
    rules <- array(sigma, c(num, 2))
    rules[,2] <- centers
  } else if (form == 'tri') {
    rules <- c(centers - center_size, centers, centers + center_size)
    rules <- matrix(rules, ncol=3)
  }
  
  mdegree <- matrix(numeric(length(y) * num), ncol=num)
  for (i in 1:num) {
    if (form == 'gauss') {
      mdegree[, i] <- exp(-(y - rules[i, 2])^2/(2*rules[i, 1]^2))
    } else if (form == 'tri') {
      mdegree[, i] <- (y > rules[i, 1] & y <= rules[i, 2])*(y - rules[i, 1])/(rules[i, 2] - rules[i, 1])
        + (y > rules[i, 2] & y <= rules[i, 3])*(rules[i, 3] - y)/(rules[i, 3] - rules[i, 2])
      
    }
  }
  
  argsortrow <- function(x) { return(sort(x, index.return=T)$ix) }  
  idx <- t(apply(mdegree, 1, argsortrow))
  return(idx[, num])
}

tslda <- function(TimeSeries, NumTopics, FuzzyForm='gauss', NumWords=5, NumWordsPerDoc=10, Iterations=300, Alpha=0.1, Gamma=0.1) {
# TSLDA  Utilização do modelo Alocação Latente de Dirichlet (LDA, do inglês
# Latent Dirichlet Allocation) para agrupamento de séries temporais.
#
# Entradas
# ========
#   Required  TimeSeries      <array <numeric>>
#     A série temporal a ser agrupada
#   Required  NumTopics       <integer>
#     O número de grupos que se espera encontrar na série
#   Optional  FuzzyForm       <string>
#     A forma fuzzy com que se deseja fuzzificar a série temporal
#     Options: gauss (default), tri
#   Optional  NumWords        <integer>
#     O número de palavras ou regras fuzzy que se deseja obter
#     Default: 5
#   Optional  NumWordsPerDoc  <integer>
#     O número de palavras que se terá em cada documento
#     Default: 10
#   Optional  Iterations      <integer>
#     O número de iterações para o amostrador de Gibbs
#     Default: 300
#   Optional  Alpha           <array <numeric>>
#     Um vetor com valores de alpha para cada tópico ou um valor único que
#     será replicado para todos os tópicos
#     Default: 0.1
#   Optional  Gamma           <array <numeric>>
#     Um vetor com valores de gamma para cada palavra ou um valor único que
#     será replicado para todos as palavras
#     Default: 0.1
#
# Saídas
# ======
#   labels <array <integer>>
#     Rótulos encontrados para cada valor da série
#   theta  <array <numeric>>
#     Thetas encontrados pelo amostrador de Gibbs em cada instante
#   beta   <array <numeric>>
#     Betas encontrados pelo amostrador de Gibbs em cada instante
#   zd     <array <integer>>
#     Zds encontrados pelo amostrador de Gibbs em cada instante

  ## Validação dos parâmetros do modelo
  validPositiveInteger <- function(x) { is.numeric(x) & x %% 1 == 0 && (x > 0) }
  validPositiveNumeric <- function(x) { is.numeric(x) & all(x > 0) }
  validFuzzyForms <- function(x) { any(c('tri', 'gauss') == x) }
  
  if (!is.numeric(TimeSeries)) {
    stop("The time series must be a numeric array")
  }
  if (!validPositiveInteger(NumTopics)) {
    stop("The number of topics must be a positive integer")
  }
  if (!validFuzzyForms(FuzzyForm)) {
    stop("The fuzzy form must be gauss or triang")
  }
  if (!validPositiveInteger(NumWords)) {
    stop("The number of words must be a positive integer")
  }
  if (!validPositiveInteger(NumWordsPerDoc)) {
    stop("The number of words per doc must be a positive integer")
  }
  if (!validPositiveInteger(Iterations)) {
    stop("The number of iterations must be a positive integer")
  }
  if (!validPositiveNumeric(Alpha)) {
    stop("The vector of alphas must be positive numeric")
  }
  if (!validPositiveNumeric(Gamma)) {
    stop("The vector of gammas must be positive numeric")
  }

  # Valida o vetor de inicial de alphas caso seja fornecido
  alpha <- Alpha
  if (length(alpha) == 1) {
    alpha <- array(Alpha, c(NumTopics))
  } else if (length(Alpha) != NumTopics) {
    stop("The alpha vector should have the same size of the amount of topics.")
  }
  
  # Valida o vetor de inicial de gammas caso seja fornecido
  gamma <- Gamma
  if (length(gamma) == 1) {
    gamma <- array(Gamma, c(NumWords))
  } else if (length(Gamma) != NumWords) {
    stop("The gamma vector should have the same size of the amount of words.")
  }
  
  ## Tranformação da série temporal
  # Fuzzifica série de entrada formando o universo de discurso
  words <- fuzzify(TimeSeries, FuzzyForm, NumWords)
  
  # Divide palavras do dicionário em documentos de tamanho igual e conta
  # quantas palavras de cada existem em cada documento
  numberOfDocs <- ceiling(length(words)/NumWordsPerDoc)
  wordCount <- array(0, c(numberOfDocs, NumWords))
  for (d in 1:numberOfDocs) {
    indexIni <- (d - 1) * NumWordsPerDoc + 1
    indexEnd <- min(indexIni + NumWordsPerDoc - 1, length(TimeSeries))
    wordsInDocD = factor(words[indexIni:indexEnd], levels=1:NumWords)
    wordCount[d,] <- table(wordsInDocD)
  }
  
  ## Inicialização dos parâmetros do modelo LDA
  # O vetor de thetas é dado pela amostra de uma Dirichlet com parametros
  # alpha
  theta <- array(0, c(NumTopics, Iterations+1))
  theta[, 1] <- dirrnd(alpha, 1)
  
  # A matriz de betas é dado por amostras de uma Dirichlet com parametros
  # gamma
  beta <- array(0, c(NumWords, NumTopics, Iterations+1))
  beta[, 1, 1] <- dirrnd(gamma, 1)
  for (i in 2:NumTopics) {
    beta[, i, 1] <- dirrnd(gamma, 1)
  }
  
  # O vetor inicial dos tópicos dos documentos é dado pela distribuição
  # multinomial cujo parâmetro é uma distribuição uniforma, inicialmente.
  zd <- array(0, c(numberOfDocs, Iterations+1))
  zd[, 1] <- sample(1:NumTopics, numberOfDocs, replace=T)
  
  ## Amostrados de Gibbs
  for (i in 2:(Iterations+1)) {
    # O novo vetor de alphas que servirá como parâmetro para a nova
    # amostragem dos thetas depende da proporção de tópicos atuais.
    topicCount <- factor(zd[, i-1], levels=1:NumTopics)
    newAlpha <- alpha + table(topicCount)
    theta[, i] <- dirrnd(newAlpha, 1)
    
    # Computa um novo gamma para cada tópico com base na ocorrência das
    # palavras em cada tópico
    for (k in 1:NumTopics) {
      newGamma <- colSums(wordCount[zd[, i-1] == k, , drop=F]) + gamma
      beta[, k, i] <- dirrnd(newGamma, 1)
    }
    
    # Reamostra os tópicos Zd com base nos vetores de thetas e betas
    for (d in 1:numberOfDocs) {
      probs <- array(0, c(NumTopics))
      for (k in 1:NumTopics) {
        probs[k] <- log(theta[k, i]) + wordCount[d, ] %*% log(beta[, k, i])
      }
      
      probs <- exp(probs)
      probs <- probs/sum(probs)
      zd[d, i] <- sample(1:NumTopics, 1, replace=T, prob=probs)
    }
  }
  
  labels <- TimeSeries * 0
  for (i in 1:length(zd[, ncol(zd)])) {
    indexIni <- (i - 1) * NumWordsPerDoc + 1
    indexEnd <- min(indexIni + NumWordsPerDoc - 1, length(TimeSeries))
    labels[indexIni:indexEnd] <- zd[i, ncol(zd)]
  }
  
  return(list(labels=labels, theta=theta, beta=beta, zd=zd))
}

## Example
y <- c(array(1, 100), array(4, 100))
y <- y + rnorm(200, 0, 1)
g <- tslda(y, 3)