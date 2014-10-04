tol = 1e-10;
[Q R E] = qr(f,0);
diagr = abs(diag(R));
rkA = find(diagr >= tol*diagr(1), 1, 'last')
n = f(:,E(1:rkA));
e = E(1:rkA);
% for i=1:size(f,2)
%     A(i) =i;
% end
% d = setdiff(A,e);