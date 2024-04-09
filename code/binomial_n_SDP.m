close all;
clear all;
m = 10;

K = 3; 

p_var = 0.1; p_level = 0.6;
p = p_level + randn(m,1)*p_var; 
p(p>0.95) = 0.95; p(p<0.05) = 0.05;

A = rand(m,m) > 0.9; A = double((A+A')/2 > 0); 
%figure; spy(A)

% form L
A = A - diag(diag(A));
D = sum(A,2);
L = diag(D) - A; disp(['true p var = ' num2str(p'*L*p)])

%form X and n
X = randn(m,K) + 5; my_beta = ones(K,1);

n = round(X*my_beta + randn(m,1));

trial = 1;

y = [];
for ii = 1:trial
y = [y, binornd(n, p)];
end
y = mean(y,2);

P = eye(m) - X*pinv(X'*X)*X';
E1 = [eye(m); zeros(m,m)];
E2 = [zeros(m,m); eye(m)];

lambda_1 = 0.2; lambda_2 = 0.2;

cvx_begin sdp
    variable Z(2*m, 2*m)  symmetric
    minimize  -2*trace(diag(y)*E2'*Z*E1)
    Z >= 0
    -E1'*Z <= 0
    -E2'*Z <= 0
    E2'*Z <= 1
cvx_end




