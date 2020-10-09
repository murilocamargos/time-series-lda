# Time-series LDA
This project uses the Latent Dirichlet Allocation (LDA) algorithm to cluster a fuzzy time-series. More details can be found on [this report](/LDA.pdf).

## LDA
LDA is a generative statistical model [[Wikipedia]](https://en.wikipedia.org/wiki/Latent_Dirichlet_allocation) proposed by David Blei, Andrew Ng and Michael I. Jordan. LDA is a three-level hierarchical Bayesian model, in which each item of a collection is modeled as a finite mixture over an underlying set of topics [[Original Paper]](https://www.researchgate.net/publication/221620547_Latent_Dirichlet_Allocation). The best example of its usage is the clustering of a text corpora into a set of topics based on their words, i.e., each document is described as a distribution of topics and each topic as a distribution of words. A very intuitive explanation about LDA can be found [here](https://towardsdatascience.com/light-on-math-machine-learning-intuitive-guide-to-latent-dirichlet-allocation-437c81220158). A very important feature of LDA is that it is able to discard topics automatically if needed, i.e., the algorithm may assign probabilities greated than zero for a smaller number of topics you originally set.

## Fuzzy Time-series
The input time-series is fuzzified with a given fuzzy form and universe of discourse. This means that each value in the time series will become a word in this universe. A sliding window with a fixed size through these words create a set of documents. Then, each document will be assigned to a set of topics with different probabilities.

# Usage
Please, use the `run.m` or `run.r` files to see the examples.