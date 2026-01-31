function G_total = egrad_U_RCG(U, H_raw, W, sim_param)
% 功能:
%   计算目标函数 J 关于接收波束矩阵 U 的复数梯度。
% 输入 (Inputs):
%   U       - [2 x K] 
%             当前的接收波束赋形矩阵
%   H_raw   - [2 x M x K]
%             原始信道张量
%   W       - [M x K] 
%             当前的发射波束赋形矩阵
%   sim_param - Struct
%             仿真参数
% 输出 (Outputs):
%   G_total - [2 x K] 
%             目标函数关于 U 的欧式梯度矩阵。
    [~, K] = size(U);
    G_total = zeros(2, K);

    mu      = sim_param.mu;
    lam2    = sim_param.lambda_sinr;
    gamma_th= sim_param.gamma_th;
    alpha   = sim_param.alpha;
    sigma2  = sim_param.sigma2;


    V_arr = zeros(2, K, K);
    for k = 1:K
        H_k = H_raw(:,:,k);          
        V_arr(:,:,k) = conj(H_k) * W; 
    end

    S_vals = zeros(K,1); I_vals = zeros(K,1); sinr_vals = zeros(K,1);

    for k = 1:K
        u_k = U(:,k);

        v_sig = V_arr(:,k,k);
        y_sig = u_k.' * v_sig;       
        S_vals(k) = abs(y_sig)^2;

        intf = 0;
        for j = 1:K
            if j ~= k
                v_int = V_arr(:,j,k);
                y_int = u_k.' * v_int; 
                intf = intf + abs(y_int)^2;
            end
        end

        I_vals(k) = intf + sigma2;
        sinr_vals(k) = S_vals(k) / I_vals(k);
    end


    exp_terms = exp(-mu * sinr_vals);
    term_main = - exp_terms / sum(exp_terms);


    term_pen = zeros(K,1);
    for k = 1:K
        Tk = gamma_th - sinr_vals(k);
        if Tk > 50/alpha
            sp = Tk; sig = 1;
        else
            val = exp(alpha*Tk);
            sp = log(1+val)/alpha;
            sig = val/(1+val);
        end
        term_pen(k) = 2 * lam2 * sp * sig * (-1);
    end

    C_k = term_main + term_pen; 


    for k = 1:K
        u_k = U(:,k);
        grad_uk = zeros(2,1);

        v_sig = V_arr(:,k,k);
        y_sig = u_k.' * v_sig;    
        coef_S = C_k(k) / I_vals(k);

        term_S = 2 * conj(v_sig) * y_sig;  
        grad_uk = grad_uk + coef_S * term_S;

        coef_I = C_k(k) * (-S_vals(k) / I_vals(k)^2);
        term_I_sum = zeros(2,1);

        for j = 1:K
            if j ~= k
                v_int = V_arr(:,j,k);
                y_int = u_k.' * v_int;     
                term_I_sum = term_I_sum + 2 * conj(v_int) * y_int;
            end
        end

        grad_uk = grad_uk + coef_I * term_I_sum;
        G_total(:,k) = grad_uk;
    end
end
