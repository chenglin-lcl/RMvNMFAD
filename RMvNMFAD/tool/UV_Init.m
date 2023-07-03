function [U, V] = UV_Init(X, k)
    rng('default');
    % 随机初始化所有的U和V
    m = size(X, 1);
    n = size(X, 2);
    U = abs(rand(m, k));
    V = abs(rand(k, n));
    [U, V] = NormalizeUV(U, V, 1, 2);
end