dirrnd <- function(thetas, n) {
  # This function samples from a Dirichlet distribution using the
  # gamma distribution sampler.
  k <- length(thetas)
  p <- matrix(numeric(n*k), nrow=n)
  for (i in 1:k) {
    p[, i] = rgamma(n, thetas[i])
  }
  return(p/matrix(rep(rowSums(p), k), nrow=n))
}

fuzzify <- function(y, form, num) {
  # This functino fuzzifies a time-series where a fixed number of
  # fuzzy rules equaly spaced are used.
  
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
  # Uses the Latent Dirichlet Allocation (LDA) to grupo a fuzzy time
  # series.
  #
  # Inputs
  # ======
  #   Required  TimeSeries      <array <numeric>>
  #     The time-series to be clustered.
  #   Required  NumTopics       <integer>
  #     The number of topics you expect to find.
  #   Optional  FuzzyForm       <string>
  #     The fuzzy form used to fuzzify the time-series.
  #     Options: gauss (default), tri
  #   Optional  NumWords        <integer>
  #     The numer of fuzzy words you wish to use.
  #     Default: 5
  #   Optional  NumWordsPerDoc  <integer>
  #     The number of words per document.
  #     Default: 10
  #   Optional  Iterations      <integer>
  #     Maximum number of iterations of the Gibbs sampler.
  #     Default: 300
  #   Optional  Alpha           <array <numeric>>
  #     A vector with alpha values for each topic or a single value that
  #     will be replicated for all topics.
  #     Default: 0.1
  #   Optional  Gamma           <array <numeric>>
  #     A vector with gammas values for each topic or a single value that
  #     will be replicated for all topics.
  #     Default: 0.1
  #
  # Outputs
  # =======
  #   labels <array <integer>>
  #     Labels found for each measurement.
  #   theta  <array <numeric>>
  #     Estimated thetas for each instant.
  #   beta   <array <numeric>>
  #     Estimated betas for each instant.
  #   zd     <array <integer>>
  #     Estimated zds for each instant.
  

  # Parameter validation
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

  # Validates the initial vector of alphas if provided
  alpha <- Alpha
  if (length(alpha) == 1) {
    alpha <- array(Alpha, c(NumTopics))
  } else if (length(Alpha) != NumTopics) {
    stop("The alpha vector should have the same size of the amount of topics.")
  }
  
  # Validates the initial vector of gammas if provided
  gamma <- Gamma
  if (length(gamma) == 1) {
    gamma <- array(Gamma, c(NumWords))
  } else if (length(Gamma) != NumWords) {
    stop("The gamma vector should have the same size of the amount of words.")
  }
  
  # Time-series transformation
  # Fuzzyfies the input series
  words <- fuzzify(TimeSeries, FuzzyForm, NumWords)
  
  # Divides the dictionary words in document of equal size and counts
  # how mny words there is in each one of them
  numberOfDocs <- ceiling(length(words)/NumWordsPerDoc)
  wordCount <- array(0, c(numberOfDocs, NumWords))
  for (d in 1:numberOfDocs) {
    indexIni <- (d - 1) * NumWordsPerDoc + 1
    indexEnd <- min(indexIni + NumWordsPerDoc - 1, length(TimeSeries))
    wordsInDocD = factor(words[indexIni:indexEnd], levels=1:NumWords)
    wordCount[d,] <- table(wordsInDocD)
  }
  
  # Initializes the LDA parameters
  # The theta vector is given by a sample of the Dirichlet distribution
  # with parameters alpha
  theta <- array(0, c(NumTopics, Iterations+1))
  theta[, 1] <- dirrnd(alpha, 1)
  
  # The beta matrix is given by a sample of the Dirichlet distribution
  # with parameters gamma
  beta <- array(0, c(NumWords, NumTopics, Iterations+1))
  beta[, 1, 1] <- dirrnd(gamma, 1)
  for (i in 2:NumTopics) {
    beta[, i, 1] <- dirrnd(gamma, 1)
  }
  
  # The initial vector of document topics is given by the multinomial
  # distribution with parameters following a uniform distribution.
  zd <- array(0, c(numberOfDocs, Iterations+1))
  zd[, 1] <- sample(1:NumTopics, numberOfDocs, replace=T)
  
  # Gibbs sampler
  for (i in 2:(Iterations+1)) {
    # The new alpha vector that will be used as parameters to sample
    # the thetas depende on the current topic proportional
    topicCount <- factor(zd[, i-1], levels=1:NumTopics)
    newAlpha <- alpha + table(topicCount)
    theta[, i] <- dirrnd(newAlpha, 1)
    
    # Compute a new gamma for each topic using the word cound on each
    # topic
    for (k in 1:NumTopics) {
      newGamma <- colSums(wordCount[zd[, i-1] == k, , drop=F]) + gamma
      beta[, k, i] <- dirrnd(newGamma, 1)
    }
    
    # Resample the topics Zd using the theta and beta vectors
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