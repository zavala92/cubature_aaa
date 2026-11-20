function [x, info] = cub_lsq_svd(A, b)
%CUB_LSQ_SVD   Least squares by SVD with scaled columns and a rank cutoff.
%   x = cub_lsq_svd(A, b) scales the columns of A to norm 1 and drops the
%   singular values below max(size(A))*eps(s(1)), the tolerance of rank.
%   [x, info] = cub_lsq_svd(A, b) also returns the retained rank, the
%   cutoff and the singular values of the scaled matrix.
nrm = sqrt(sum(abs(A).^2, 1));
[U, S, V] = svd(A./nrm, 'econ');
s = diag(S);
tol = max(size(A))*eps(s(1));
k = s > tol;
x = (V(:,k)*((U(:,k)'*b)./s(k)))./nrm.';

info.rank = nnz(k);
info.ncols = size(A, 2);
info.tol = tol;
info.sigma = s;
info.condEff = s(1)/s(info.rank);
end
