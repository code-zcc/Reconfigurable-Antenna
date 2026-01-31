function G_total = egrad_RCG_Multipath(R, H_data, W, param)
% 目的:
%   计算目标函数 J 关于旋转矩阵 R 的欧式梯度。
% 输入:
%   R      - [3 x 3 x M] 当前的天线旋转矩阵
%   H_data - Struct      预计算的信道几何数据 (方向 d, 极化 psi 等)
%   W      - [M x K]     当前的发射波束
%   param  - Struct      仿真参数
%
% 输出:
%   G_total - [3 x 3 x M] 梯度张量
    [~, K] = size(W);
    M = size(R, 3);
    G_total = zeros(3, 3, M);
    

    [H_val, G_vals_all, P_vals_all] = compute_H_internal(R, H_data, param);
    

    S_vals = zeros(K,1); I_vals = zeros(K,1); sinr_vals = zeros(K,1);
    for k=1:K
        S_vals(k) = abs(W(:,k)' * H_val(:,k))^2;
        intf = 0;
        for j=1:K, if j~=k, intf = intf + abs(W(:,j)' * H_val(:,k))^2; end, end
        I_vals(k) = intf + param.sigma2;
        sinr_vals(k) = S_vals(k) / I_vals(k);
    end
    

    mu = param.mu;
    exp_terms = exp(-mu * sinr_vals);
    sum_exp = sum(exp_terms);

    term_main = - exp_terms / sum_exp;
    

    lam_sinr = param.lambda_sinr;
    gamma_th = param.gamma_th;
    alpha = param.alpha;
    term_pen = zeros(K,1);
    for k=1:K
        delta = gamma_th - sinr_vals(k);
        if delta > 50/alpha
            sp = delta; sig = 1;
        else
            val = exp(alpha*delta);
            sp = (1/alpha)*log(1+val);
            sig = val/(1+val);
        end

        term_pen(k) = 2 * lam_sinr * sp * sig * (-1);
    end
    C_k = term_main + term_pen;
    

    z_hat = [0;0;1];
    cos_max = cos(param.theta_max);
    lam_ang = param.lambda_angle;
    
    for m=1:M
        grad_Rm = zeros(3,3);
        r1=R(:,1,m); r2=R(:,2,m); r3=R(:,3,m);
        
        for k=1:K

            dh_dr1 = zeros(3,1); dh_dr2 = zeros(3,1); dh_dr3 = zeros(3,1);
            
            for l=1:H_data.num_paths
                d_vec = H_data.d_vecs(:,k,m,l);
                psi   = H_data.psi_vecs(:,k,m,l);
                C0    = H_data.consts(k,m,l);
                G     = G_vals_all(m,k,l);
                P     = P_vals_all(m,k,l);
                p     = param.p;
                fv    = H_data.f_analog(1,m);
                fh    = H_data.f_analog(2,m);
                

                dh_dr1 = dh_dr1 + C0 * G * (fv * conj(psi));
                dh_dr2 = dh_dr2 + C0 * G * (fh * conj(psi));
                
                u_dot = d_vec' * r3;
                if u_dot > 0
                    dh_dr3 = dh_dr3 + C0 * P * p * (u_dot)^(p-1) * d_vec;
                end
            end
            
            coef_S = C_k(k) / I_vals(k);
            y_sig = W(:,k)' * H_val(:,k);
            scalar_S = conj(y_sig) * conj(W(m,k));
            
            grad_Rm(:,1) = grad_Rm(:,1) + coef_S * 2*real(scalar_S * dh_dr1);
            grad_Rm(:,2) = grad_Rm(:,2) + coef_S * 2*real(scalar_S * dh_dr2);
            grad_Rm(:,3) = grad_Rm(:,3) + coef_S * 2*real(scalar_S * dh_dr3);
            

            coef_I = C_k(k) * (-S_vals(k) / I_vals(k)^2);
            for j=1:K
                if j~=k
                    y_intf = W(:,j)' * H_val(:,k);
                    scalar_I = conj(y_intf) * conj(W(m,j));
                    
                    grad_Rm(:,1) = grad_Rm(:,1) + coef_I * 2*real(scalar_I * dh_dr1);
                    grad_Rm(:,2) = grad_Rm(:,2) + coef_I * 2*real(scalar_I * dh_dr2);
                    grad_Rm(:,3) = grad_Rm(:,3) + coef_I * 2*real(scalar_I * dh_dr3);
                end
            end
        end
        

        delta = cos_max - r3'*z_hat;

             if delta > 50/alpha
                 sp = delta; sig = 1;
             else
                 val = exp(alpha*delta);
                 sp = (1/alpha)*log(1+val);
                 sig = val/(1+val);
             end
             grad_ang = 2 * lam_ang * sp * sig * (-z_hat);
             grad_Rm(:,3) = grad_Rm(:,3) + grad_ang;

        
        G_total(:,:,m) = grad_Rm;
    end

end