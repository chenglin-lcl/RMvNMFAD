% compute squared Euclidean distance
% ||A-B||^2 = ||A||^2 + ||B||^2 - 2*A'*B
function d = L2_distance_1(a,b)
% a,b: two matrices. each column is a data
% d:   distance matrix of a and b
% dij是指a的第i列与b的第j列之间的欧式距离的平方


if (size(a,1) == 1) % 如果a的行数是1
  a = [a; zeros(1,size(a,2))]; % 在a的下面补个0行
  b = [b; zeros(1,size(b,2))]; % 在a的下面补个0行
end

aa=sum(a.*a); % a中的所有元素平方
bb=sum(b.*b); 
ab=a'*b;  
d = repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab;

d = real(d);
d = max(d,0);


% % force 0 on the diagonal? 
% if (df==1)
%   d = d.*(1-eye(size(d)));
% end
