function G_r3 = egrad_RCG_PointingOnly(r3_list, r1_fixed, r2_fixed, H_data, W, param)
    % 输入 R: [3 x 3 x M] 组装好的矩阵 (包含固定的 r1, r2 和 变量 r3)
    % 输出 G_r3: [3 x M] 仅包含对 r3 的欧氏梯度
    [~, M] = size(r3_list);
    
    % =========================================================
    % 1. 直接拼接 (Direct Assembly)
    % =========================================================
    % 此时 R 可能暂时不满足刚体正交约束，但这正是你AO策略里“松弛”的一步
    R = zeros(3, 3, M);
    for m = 1:M
        R(:, 1, m) = r1_fixed(:, m); % 锁死
        R(:, 2, m) = r2_fixed(:, m); % 锁死
        R(:, 3, m) = r3_list(:, m);  % 变量
    end
    [~, K] = size(W);
%     M = size(R, 3);
    G_r3 = zeros(3, M); % 输出维度变了
    
    % Re-calculate Forward Pass (需要完整的 R 来算准确的 SINR)
    [H_val, G_vals_all, P_vals_all] = compute_H_internal(R, H_data, param);
    
    % --- 1. Calculate SINR & Common Terms (不变) ---
    S_vals = zeros(K,1); I_vals = zeros(K,1); sinr_vals = zeros(K,1);
    for k=1:K
        S_vals(k) = abs(W(:,k)' * H_val(:,k))^2;
        intf = 0;
        for j=1:K, if j~=k, intf = intf + abs(W(:,j)' * H_val(:,k))^2; end, end
        I_vals(k) = intf + param.sigma2;
        sinr_vals(k) = S_vals(k) / I_vals(k);
    end
    
    % Coefficient C_k (Chain Rule part)
    mu = param.mu;
    exp_terms = exp(-mu * sinr_vals);
    sum_exp = sum(exp_terms);
    term_main = - exp_terms / sum_exp; % Soft-min derivation
    
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
    
    % --- 2. Backpropagation (只算 r3) ---
    z_hat = [0;0;1];
    cos_max = cos(param.theta_max);
    lam_ang = param.lambda_angle;
    
    for m=1:M
        grad_r3_m = zeros(3,1); % 单个天线 r3 的梯度
        r3 = R(:,3,m);          % 我们只需要 r3
        
        for k=1:K
            % === 关键修改：只计算 dh/dr3 ===
            % dh/dr1 和 dh/dr2 直接删掉，不算了，省时间
            dh_dr3 = zeros(3,1);
            
            for l=1:H_data.num_paths
                d_vec = H_data.d_vecs(:,k,m,l);
                % psi, fv, fh 这些跟 r1, r2 有关的项，对于 r3 的导数没贡献，跳过
                
                % r3 只影响增益项 (Gain Term)
                % 注意：如果你的模型里 P (Polarization term) 也显式包含了 r3，
                % 那么 P 对 r3 的导数也要算。
                % 但通常 P = u' * Q * R * a，如果 a=[1;0;0] (本地坐标)，
                % 那么 P 其实只跟 r1, r2 有关。
                % **这里假设 r3 主要影响 Directivity Gain (G)**
                
                C0    = H_data.consts(k,m,l);
                P     = P_vals_all(m,k,l); % P 值还要用到，因为它在链式法则系数里
                p     = param.p;
                
                u_dot = d_vec' * r3;
                if u_dot > 0
                    % 这是 Gain 对 r3 的导数
                    dh_dr3 = dh_dr3 + C0 * P * p * (u_dot)^(p-1) * d_vec;
                end
            end
            
            % Map to J via SINR
            % Signal Part
            coef_S = C_k(k) / I_vals(k);
            y_sig = W(:,k)' * H_val(:,k);
            scalar_S = conj(y_sig) * conj(W(m,k));
            
            grad_r3_m = grad_r3_m + coef_S * 2*real(scalar_S * dh_dr3);
            
            % Interference Part
            coef_I = C_k(k) * (-S_vals(k) / I_vals(k)^2);
            for j=1:K
                if j~=k
                    y_intf = W(:,j)' * H_val(:,k);
                    scalar_I = conj(y_intf) * conj(W(m,j));
                    grad_r3_m = grad_r3_m + coef_I * 2*real(scalar_I * dh_dr3);
                end
            end
        end
        
        % 3. Angle Penalty (只针对 r3，必须保留)
        delta = cos_max - r3'*z_hat;
        if delta > 50/alpha
            sp = delta; sig = 1;
        else
            val = exp(alpha*delta);
            sp = (1/alpha)*log(1+val);
            sig = val/(1+val);
        end
        grad_ang = 2 * lam_ang * sp * sig * (-z_hat);
        
        % 累加 Penalty 梯度
        grad_r3_m = grad_r3_m + grad_ang;
        
        % 存入结果
        G_r3(:, m) = grad_r3_m;
    end
end