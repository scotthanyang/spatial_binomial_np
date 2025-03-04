

close all;
clear all;

rng(2024);
%%
m = 10;

K = 3; 

p_var = 0.3; p_level = 0.5;
p = p_level + randn(m,1)*p_var; 
p(p>0.95) = 0.95; p(p<0.05) = 0.05;

% make p bounded between [p_min, p_max]
p_min = 0.05;
p_max = 0.95;
p = p_min+p*(p_max-p_min);

A = rand(m,m) > 0.9; A = double((A+A')/2 > 0);

A = A - diag(diag(A));
D = sum(A,2);
L = diag(D) - A; 
disp(['true p var = ' num2str(p'*L*p)])

%%
X = randn(m,K) + 5; my_beta = ones(K,1);
n = round(X*my_beta + randn(m,1));
%X = randn(m,K) + 2; my_beta = ones(K,1);
%n = round(exp(X*my_beta+ randn(m,1)));


trial = 10;
tic,
y = [];
for ii = 1:trial
    y = [y, binornd(n, p)];
end
y = mean(y,2);
toc

%%
%lambda_1 = 0.2; 
lambda_1 = 0.01; 

lambda_2 = 0.5;
P = eye(m) - X*pinv(X'*X)*X';

%%
iter = 3000; %500; 
loss = zeros(iter,1);
loss2 = zeros(iter,1); %log domain in 1st term


% initial value
n_est = y + 1; 

%p_est = ones(m,1)*p_level;
p_est = y./n_est;
%p_est = p_min + rand(m,1)*(p_max-p_min);


%% in (log n, log p) domain

v_min = log(p_min);
v_max = log(p_max);

u = log(n_est);
v = log(p_est);

logy = log(y);


% iterations
step = 1e-5; %1e-3; 
inner_iter = 1;

for ii = 1:iter
    
    for inter_ii = 1:inner_iter
        
        %dn1 = -2*(y-n_est.*p_est).*p_est ;
        du1 = -2* (logy -(u + v));
        
        n_est = exp(u);
        dn2 = 2*lambda_2*P*n_est;
        du2 = dn2.*n_est;
        
        u = u - step*(du1+du2);
        u_min = 0;
        
        u( u < u_min) =u_min;
        
        %n_est = n_est - step*(dn1+dn2);
        %n_est = n_est - step*(-2*(y-n_est.*p_est).*p_est + 2*lambda_2*P*n_est);
        
        
        %n_est(n_est<1) =1;
    end
    
    for inter_ii = 1:inner_iter
        
        %dp1 = -2*(y-n_est.*p_est).*n_est;
        dv1 = -2*(logy -(u+v));
        
        p_est = exp(v);
        dp2 =  2*lambda_1*L*p_est;
        dv2 = dp2.*p_est;
        
        v = v - step*(dv1+dv2);
        
        v(v>v_max) = v_max; 
        v(v<v_min) = v_min;
        
        %p_est = p_est - step*( dp1 + dp2 );
        
        %p_est = p_est - step*( -2*(y-n_est.*p_est).*n_est + 2*lambda_1*L*p_est );
        %p_est(p_est>1) = 1; p_est(p_est<0) = 0;
        
        %p_est(p_est>p_max) = p_max; 
        %p_est(p_est<p_min) = p_min;
    end
    
    n_est = exp(u);
    p_est = exp(v);
    
    loss(ii) = norm(y - n_est.*p_est).^2 + lambda_1*p'*L*p ...
        + lambda_2*norm(n-X*my_beta).^2;
    
    loss2(ii) = norm(logy - (u+v)).^2 + lambda_1*p'*L*p ...
        + lambda_2*norm(n-X*my_beta).^2;
    
end


%%
figure(1);clf; 
subplot(121);
plot(loss,'-o'); 
title('Loss'); grid on;
subplot(122);
plot(loss2,'-o'); 
title('Loss2'); grid on;

figure; plot(loss,'-*'); ylabel('Loss'); xlabel('iteration')

figure; subplot(2,1,1); stem(p); hold on; stem(p_est,'r'); legend('true p','est p')

subplot(2,1,2); stem(n); hold on; stem(n_est,'r');  stem(y,'g'); legend('true n','est n','y')

disp('n, y, n_est')
[n y n_est]
% %%
% disp(' p, p_est')
% [p, p_est]
% 
% disp(' n, n_est')
% [n, n_est]

%return;
%


% %% in (n,p) domain
% 
% % iterations
% step = 0.01; 
% inner_iter = 7;
% 
% for ii = 1:iter
%     
%     for inter_ii = 1:inner_iter
%         dn1 = -2*(y-n_est.*p_est).*p_est ;
%         dn2 = 2*lambda_2*P*n_est;
%         n_est = n_est - step*(dn1+dn2);
%         %n_est = n_est - step*(-2*(y-n_est.*p_est).*p_est + 2*lambda_2*P*n_est);
%         
%         n_est(n_est<1) =1;
%     end
%     
%     for inter_ii = 1:inner_iter
%         
%         dp1 = -2*(y-n_est.*p_est).*n_est;
%         dp2 =  2*lambda_1*L*p_est;
%         p_est = p_est - step*( dp1 + dp2 );
%         
%         %p_est = p_est - step*( -2*(y-n_est.*p_est).*n_est + 2*lambda_1*L*p_est );
%         %p_est(p_est>1) = 1; p_est(p_est<0) = 0;
%         
%         p_est(p_est>p_max) = p_max; 
%         p_est(p_est<p_min) = p_min;
%     end
%     
%     loss(ii) = norm(y - n_est.*p_est).^2 + lambda_1*p_est'*L*p_est ...
%         + lambda_2*norm(n_est-X*my_beta).^2;
%     %loss(ii)
% end
% 
% %%
% figure; plot(loss,'-o'); title('Loss')
% 
% figure; subplot(2,1,1); stem(p); hold on; stem(p_est,'r'); legend('true p','est p')
% 
% subplot(2,1,2); stem(n); hold on; stem(n_est,'r');  stem(y,'g'); legend('true n','est n','y')
% 
% disp('n, y, n_est')
% [n y n_est]
