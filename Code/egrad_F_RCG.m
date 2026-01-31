function G_total = egrad_F_RCG(F, R_curr, H_func_data, W, sim_param)
    [~, K] = size(W);
    [~, M] = size(F);
    G_total = zeros(2, M);
    

    mu = sim_param.mu;
    lam2 = sim_param.lambda_sinr;
    gamma_th = sim_param.gamma_th;
    alpha = sim_param.alpha;
    

    [H_val, G_vals_all] = compute_H_with_F(R_curr, F, H_func_data, sim_param);
    

    S_vals = zeros(K,1); I_vals = zeros(K,1); sinr_vals = zeros(K,1);
    for k=1:K
        S_vals(k) = abs(W(:,k)' * H_val(:,k))^2;
        intf = 0;
        for j=1:K, if j~=k, intf = intf + abs(W(:,j)' * H_val(:,k))^2; end, end
        I_vals(k) = intf + sim_param.sigma2;
        sinr_vals(k) = S_vals(k) / I_vals(k);
    end
    

    exp_terms = exp(-mu * sinr_vals);
    term_main = - exp_terms / sum(exp_terms);
    
    term_pen = zeros(K,1);
    for k=1:K
        Tk = gamma_th - sinr_vals(k);
        if Tk > 50/alpha, sp=Tk; sig=1; else, val=exp(alpha*Tk); sp=log(1+val)/alpha; sig=val/(1+val); end
        term_pen(k) = 2 * lam2 * sp * sig * (-1);
    end
    C_k = term_main + term_pen;
    

    for m = 1:M
        grad_fm = zeros(2, 1);
        r1 = R_curr(:,1,m);
        r2 = R_curr(:,2,m);
        
        for k = 1:K
            
            beta_vec = zeros(1, 2);
            for l = 1:H_func_data.num_paths
                psi_stored = H_func_data.psi_vecs(:, k, m, l);
                C0 = H_func_data.consts(k, m, l);
                G  = G_vals_all(m, k, l);
                common = C0 * G;
               
                beta_vec(1) = beta_vec(1) + common * (r1'*conj(psi_stored));
                beta_vec(2) = beta_vec(2) + common * (r2'*conj(psi_stored));
            end
            

            dh_df_conj = beta_vec'; % 2x1 Column Vector

            coef_S = C_k(k) / I_vals(k);
            y_sig = W(:,k)' * H_val(:,k);%H_val(:,k)'*W(:,k);
            scalar_S = (y_sig) * W(m,k); %conj(y_sig) * W(m,k); 
            
            grad_fm = grad_fm + coef_S * 2 * scalar_S * dh_df_conj;
            
            coef_I = C_k(k) * (-S_vals(k) / I_vals(k)^2);
            term_I = zeros(2,1);
            for j = 1:K
                if j~=k
                    y_intf =  W(:,j)'*H_val(:,k);%H_val(:,k)'*W(:,j);
                    scalar_I = (y_intf) *(W(m,j));
                    term_I = term_I + 2 * scalar_I * dh_df_conj;
                end
            end
            grad_fm = grad_fm + coef_I * term_I;
        end
        G_total(:, m) = grad_fm;
    end
end