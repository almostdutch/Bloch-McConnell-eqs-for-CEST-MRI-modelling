function value = matrix_exp(A)
% matrix exponential

[V,D] = eig(A);
value = V*diag(exp(diag(D)))*V';
end