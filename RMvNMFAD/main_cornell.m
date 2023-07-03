clear;
clc;

addpath('tool/', 'update/', 'dataset/');

load cornell.mat
X{1} = X1; 
X{2} = X2; 
truth = Y'; 

option = [];
option.view_num = length(X); % 视图个数
option.maxIter = option.view_num * 500;
option.k = length(unique(truth)); % 簇的个数 ||X-UV|| U(m*k) V(k*n)

option.alpha = 1000;
option.gamma = 0.001;
option.lambda = 0.1;

[V_star, obj] = MY_MVC_update(X, option);
[AC, NMI, purity] = calc_Metric(V_star', option.k, truth, 10);

plot(obj, 'LineWidth', 1, 'Color',[0 0 1]);
xlabel('Number of iteration');
ylabel('Objective function value');
title('Cornell');




