%% UV单位化
function [U, V] = NormalizeUV(U, V, NormV, Norm)
% 适用于 X-UV
    K = size(U,2);
    if Norm == 2
        if NormV
            norms = max(1e-15,sqrt(sum(V.^2,1)))';
            V = spdiags(norms.^-1,0,K,K)*V;
            U = U*spdiags(norms,0,K,K);
        else 
            % 对U m*k进行列分块，U = {u1, u2, ..., uk}
            % Q = diag{||u1||, ||u2||, ..., ||uk||}
            norms = max(1e-15,sqrt(sum(U.^2,1)))';
            % U = U*inv(Q) = {u1/||u1||, u2/||u2||, ..., uk/||uk||}
            U = U*spdiags(norms.^-1,0,K,K);
            V = spdiags(norms,0,K,K)*V;
            
            %X = UV => X = UQ^(-1)QV
        end
    else
        if NormV
            norms = max(1e-15,sum(abs(V),1))';
            V = spdiags(norms.^-1,0,K,K)*V;
            U = U*spdiags(norms,0,K,K);
        else
            norms = max(1e-15,sum(abs(U),1))';
            U = U*spdiags(norms.^-1,0,K,K);
            V = spdiags(norms,0,K,K)*V;
        end
    end
end
