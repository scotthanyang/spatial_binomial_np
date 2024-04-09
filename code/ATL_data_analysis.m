

close all;
clear all;

%load atl_data.mat;

load atl_data_no88.mat; y = y/500;

%%
m = 78;

K = 3; 

% make p bounded between [p_min, p_max]
p_min = 0.1;
p_max = 0.9;

A = A - diag(diag(A));
D = sum(A,2);
L = diag(D) - A; 


%%
%lambda_1 = 0.2; 
lambda_1 = 0.001; 
lambda_2 = 0.001;

X = [ones(m,1), X]; K = K+1;
P = eye(m) - X*pinv(X'*X)*X';

%%
iter = 100000; %500; 
loss = zeros(iter,1);
loss2 = zeros(iter,1); %log domain in 1st term


% initial value
%n_est = y + 20;  p_est = y./n_est;
p_level = 0.8; p_est = ones(m,1)*p_level; n_est = y./p_est;

%n_est = y./p_est;
%p_est = y./n_est;
%p_est = p_min + rand(m,1)*(p_max-p_min);


%     %% in (log n, log p) domain
%     
%     v_min = log(p_min);
%     v_max = log(p_max);
%     
%     u = log(n_est);
%     v = log(p_est);
%     
%     logy = log(y);
%     
%     
%     % iterations
%     step = 1e-2; %1e-3; 
%     inner_iter = 5;
%     y_scale = 1/m;
%     
%     for ii = 1:iter
%         
%         for inter_ii = 1:inner_iter
%             
%             %dn1 = -2*(y-n_est.*p_est).*p_est ;
%             du1 = -2* (logy -(u + v))*y_scale;
%             
%             n_est = exp(u);
%             dn2 = 2*lambda_2*P*n_est;
%             du2 = dn2.*n_est;
%             
%             u = u - step*(du1+du2);
%             u_min = 0;
%             
%             u( u < u_min) =u_min;
%             
%             %n_est = n_est - step*(dn1+dn2);
%             %n_est = n_est - step*(-2*(y-n_est.*p_est).*p_est + 2*lambda_2*P*n_est);
%             
%             
%             %n_est(n_est<1) =1;
%         end
%         
%         for inter_ii = 1:inner_iter
%             
%             %dp1 = -2*(y-n_est.*p_est).*n_est;
%             dv1 = -2*(logy -(u+v))*y_scale;
%             
%             p_est = exp(v);
%             dp2 =  2*lambda_1*L*p_est;
%             dv2 = dp2.*p_est;
%             
%             v = v - step*(dv1+dv2);
%             
%             v(v>v_max) = v_max; 
%             v(v<v_min) = v_min;
%             
%             %p_est = p_est - step*( dp1 + dp2 );
%             
%             %p_est = p_est - step*( -2*(y-n_est.*p_est).*n_est + 2*lambda_1*L*p_est );
%             %p_est(p_est>1) = 1; p_est(p_est<0) = 0;
%             
%             %p_est(p_est>p_max) = p_max; 
%             %p_est(p_est<p_min) = p_min;
%         end
%         
%         n_est = exp(u);
%         p_est = exp(v);
%         
%         loss(ii) = (1/m.^2)*norm(y - n_est.*p_est).^2 + lambda_1*p_est'*L*p_est ...
%             + lambda_2*n_est'*P*n_est;
%         
%         %loss2(ii) = norm(logy - (u+v)).^2 + lambda_1*p'*L*p ...
%             %+ lambda_2*norm(n-X*my_beta).^2;
%         
%     end
%     
%     
%     %%
%     figure(1);%clf; 
%     %subplot(121);
%     plot(loss,'-o'); 
%     title('Loss'); grid on;
%     % subplot(122);
%     % plot(loss2,'-o'); 
%     % title('Loss2'); grid on;
%     
%     figure; subplot(2,1,1);  stem(p_est,'r'); legend('est p')
%     
%     subplot(2,1,2); stem(y); hold on; stem(n_est,'r'); legend('y','est n')
%     
%     disp('n, y, n_est')
%     [y n_est]
%     % %%
%     % disp(' p, p_est')
%     % [p, p_est]
%     % 
%     % disp(' n, n_est')
%     % [n, n_est]
%     
%     return;
%     %%


%% in (n,p) domain

% iterations
step = 0.001; 
inner_iter = 10;
p_diff=[];
n_diff=[];

for ii = 1:iter
    
    for inter_ii = 1:inner_iter
        dn1 = -2*(y-n_est.*p_est).*p_est/m;
        dn2 = 2*lambda_2*P*n_est;
        n_est = n_est - step*(dn1+dn2);
        %n_est = n_est - step*(-2*(y-n_est.*p_est).*p_est + 2*lambda_2*P*n_est);
        
        n_est(n_est<1) =1;
    end
    
    for inter_ii = 1:inner_iter
        
        dp1 = -2*(y-n_est.*p_est).*n_est/m;
        dp2 =  2*lambda_1*L*p_est;
        p_est = p_est - step*( dp1 + dp2 );
        
        %p_est = p_est - step*( -2*(y-n_est.*p_est).*n_est + 2*lambda_1*L*p_est );
        %p_est(p_est>1) = 1; p_est(p_est<0) = 0;
        
        p_est(p_est>p_max) = p_max; 
        p_est(p_est<p_min) = p_min;
    end
    
    loss(ii) = norm(y - n_est.*p_est).^2/m + lambda_1*p_est'*L*p_est ...
        + lambda_2*n_est'*P*n_est;
    %loss(ii)
    if ii > 1
        p_diff = [p_diff norm(dp1 + dp2)];
        n_diff = [n_diff norm(dn1 + dn2)];
    end

end

%%
figure; plot(loss,'-'); title('Loss')

figure; subplot(2,1,1); stem(p_est); legend('est p')

subplot(2,1,2); stem(n_est); hold on; stem(y,'g'); legend('est n','y')

figure; subplot(2,1,1); plot(p_diff);
subplot(2,1,2); plot(n_diff)

disp('n, y, n_est')
[y n_est]

%%loc

[val, loc] = sort(p_est);

beat(loc(1:5))

figure; stem(sort(p_est))







