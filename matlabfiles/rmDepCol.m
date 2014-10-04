A=[];
d=[];
%in=[1 2 3;4 5 6;7 8 9];
tol = 1e-10;
[Q, R, E] = qr(in);
diagr = abs(diag(R));
rkA = find(diagr >= tol*diagr(1), 1, 'last')
e = E(1:rkA);
for i=1:numPBs
    A(i) =i;
end
d = setdiff(A,e);
