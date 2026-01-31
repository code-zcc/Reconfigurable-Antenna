function [R_opt, min_sinr_val] = run_RCG_step_pointing_only(r3_init, r1_fixed, r2_fixed, W, H_func_data, sim_param)
    % 输入:
    %   r3_init:  [3 x M] 初始值 (只包含第三列)
    %   r1_fixed: [3 x M] 固定的第一列 (常量)
    %   r2_fixed: [3 x M] 固定的第二列 (常量)
    %   W, H_func_data, sim_param: 其他参数
    
    % 输出:
    %   R_opt: [3 x 3 x M] 优化后的完整矩阵 (拼好的)


    problem.M = obliquefactory(3, sim_param.M);
    

    problem.cost = @(r3_curr) calc_objective_fixed_r12(r3_curr, r1_fixed, r2_fixed, ...
                                                       H_func_data, W, sim_param);
    

    problem.egrad = @(r3_curr) egrad_RCG_PointingOnly(r3_curr, r1_fixed, r2_fixed, H_func_data, W, sim_param);
                                                 

    options.maxiter = 50; 
    options.tolgradnorm = 1e-3;
    options.verbosity = 0; 
    options.beta_type = 'P-R';
    options.linesearch = @linesearch_adaptive;
    grad_euclidean = problem.egrad(r3_init);
    grad_riemannian = problem.M.egrad2rgrad(r3_init, grad_euclidean);
  
    r3_opt = conjugategradient(problem, r3_init, options);
    

    R_opt=rotate_frame_z_to_r3(r3_opt);
    r1_fixed = R_opt(:, 1, :); 
        r2_fixed = R_opt(:, 2, :);
        r3_curr  = R_opt(:, 3, :);
        r3_curr  = squeeze(r3_curr);
    r1_fixed = squeeze(r1_fixed);
    r2_fixed = squeeze(r2_fixed);

    [~, sinr_vals] = calc_objective_fixed_r12(r3_curr, r1_fixed, r2_fixed, ...
                                              H_func_data, W, sim_param);
    min_sinr_val = min(sinr_vals);
end