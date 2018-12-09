% Série temporal
y = [ones(1, 100) ones(1, 100) * 5];
y = y + randn(1, 200);

%% LDA
[labels, theta, beta, zd] = tslda(y, 3);
unique(labels)