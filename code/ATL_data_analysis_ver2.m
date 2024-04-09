

close all;
clear all;

load atl_data.mat;

%%
m = 78;

K = 3; 


% make p bounded between [p_min, p_max]
p_min = 0.01;
p_max = 0.99;

A = A - diag(diag(A));
D = sum(A,2);
L = diag(D) - A; 


%%
%lambda_1 = 0.2; 
lambda_1 = 0.001; 
lambda_2 = 0.001;

P = eye(m) - X*pinv(X'*X)*X';

%%
iter = 250; %500; 
loss = zeros(iter,1);
loss2 = zeros(iter,1); %log domain in 1st term


% initial value


p_level = 0.5;
p_est = ones(m,1)*p_level;
%p_est = y./n_est;
%p_est = p_min + rand(m,1)*(p_max-p_min);

n_est = y./p_est; 


    % in (log n, log p) domain
    
    v_min = log(p_min);
    v_max = log(p_max);
    
    u = log(n_est);
    v = log(p_est);
    
    logy = log(y);
    
    
    %iterations
    step = 1e-2; %1e-3; 
    inner_iter = 5;
    y_scale = 1/m;
    
    for ii = 1:iter
        
        for inter_ii = 1:inner_iter
           
            du1 = -2* (logy -(u + v))*y_scale;
            
            n_est = exp(u);
            dn2 = 2*lambda_2*P*n_est;
            du2 = dn2.*n_est;
            
            u = u - step*(du1+du2);
            u_min = 0;
            
            u(u < u_min) =u_min;
            
        end
        
        for inter_ii = 1:inner_iter
            
            dv1 = -2*(logy -(u+v))*y_scale;
            
            p_est = exp(v);
            dp2 =  2*lambda_1*L*p_est;
            dv2 = dp2.*p_est;
            
            v = v - step*(dv1+dv2);
            
            v(v>v_max) = v_max; 
            v(v<v_min) = v_min;
            
        end
        
        n_est = exp(u);
        p_est = exp(v);
        
        loss(ii) = y_scale*norm(y - n_est.*p_est).^2 + lambda_1*p_est'*L*p_est ...
            + lambda_2*n_est'*P*n_est;
        
        loss2(ii) = norm(logy - (u+v)).^2 + lambda_1*p_est'*L*p_est ...
            + lambda_2*n_est'*P*n_est;
        
    end
    
    n_est = y./p_est;
    %
    figure(1);%clf; 
    subplot(121);
    plot(loss,'-o'); 
    title('Loss'); grid on;
    subplot(122);
    plot(loss2,'-o'); 
    title('Loss2'); grid on;
    
    figure; subplot(2,1,1);  stem(p_est,'r'); legend('est p')
    
    subplot(2,1,2); stem(y); hold on; stem(n_est,'r'); legend('y','est n')
    
    disp('y, n_est')
    [y n_est]
    %%
    disp('p, p_est')
    [p_est]








