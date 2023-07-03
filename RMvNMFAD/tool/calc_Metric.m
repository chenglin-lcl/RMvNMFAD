% Run k-means n times and report means and standard deviations of the
% performance measures.
%
% -------------------------------------------------------
% Input:
%       X:  data matrix (rows are samples) 行是样本
%       k:  number of clusters
%       truth:  truth cluster indicators 行向量
%     
%
% Output:
%       AC:  clustering accuracy (mean +stdev)
%       nmi:  normalized mutual information (mean +stdev)
%       purity
%
function [AC, NMI, purity] = calc_Metric(X, k, truth, replic)

if (min(truth)==0)
    truth = truth+1;
end
rng('default');
AC_ = zeros(1, replic);
NMI_ = zeros(1, replic);
purity_ = zeros(1, replic);
for i=1:replic
    idx = litekmeans(X, k, 'Replicates', 20);
    result = ClusteringMeasure(truth, idx);
    AC_(i) = result(1);
    NMI_(i) = result(2);
    purity_(i) = result(3);
end
% 求每个指标均值和方差
AC(1) = mean(AC_); AC(2) = std(AC_);
NMI(1) = mean(NMI_); NMI(2) = std(NMI_);
purity(1) = mean(purity_); purity(2) = std(purity_);

% 打印结果
fprintf("AC = %5.4f + %5.4f, NMI = %5.4f + %5.4f, purity = %5.4f + %5.4f\n", AC(1), AC(2), NMI(1), NMI(2), purity(1), purity(2));

end