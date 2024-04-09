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

% 
% Z = randn(2*m, 2*m);
% Z = (Z+Z')/2;
% [U, D] = eig(Z);
% Z = U*diag(abs(diag(D)))*U';

iter = 10; loss = zeros(iter,1);
inner_iter = 10;

p_level = 0.8; p_est = p_level*ones(m,1);
n_est = y./p_est;
z = [n_est; p_est];
Z = z*z';

EE = E2*E1';
step = 0.01;
for ii = 1:iter
    
    d1 = 2*EE'*Z*EE';
    d2 = -2*E2*diag(y)*E1';
    d3 = E2*L*E2';
    d4 = E1*P*E1';

    Z = Z - step*(d1 + d2 + lamda_1*d3 + lambda_2*d4);

    [U, D] = eig(Z); d = diag(D); d = d.*(d>=0.01);

    loss(ii) = trace(E1'*Z*E2*E1'*Z*E2) - 2;
end
