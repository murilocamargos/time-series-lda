clear; clc;
close all;

addpath('src');

%% Create a time series with an abrupt changepoint
y = [ones(1, 100) ones(1, 100) * 5];
y = y + randn(1, 200);

%% Run TSLDA with an approximated number of clusters (3)
[labels, theta, beta, zd] = tslda(y, 3);