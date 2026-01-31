function [R_opt, min_sinr_val] = run_RCG_step(R_init, W, H_func_data, sim_param)
    % 1. 准备 Manopt 问题结构
    % 我们的变量是 M 个 3x3 旋转矩阵，属于 SO(3)^M 流形
    problem.M = rotationsfactory(3, sim_param.M);
    
    % 2. 定义目标函数 (Cost)
    % 注意：H_func_data 包含了冻结的信道几何信息
    problem.cost = @(R) calc_objective(R, H_func_data, W, sim_param);
    
    % 3. 定义欧式梯度 (Euclidean Gradient)
    % Manopt 会自动将其转换为黎曼梯度
    problem.egrad = @(R) egrad_RCG_Multipath(R, H_func_data, W, sim_param);
    
    % 4. 配置求解器选项
    options.maxiter = 50;        % 内部迭代次数，不需要太多，因为外面还有 AO 循环
    options.tolgradnorm = 1e-3;  % 梯度收敛门限
    options.verbosity = 0;       % 0: 静默模式 (不刷屏), 1: 显示日志
    options.beta_type = 'P-R';   % Polak-Ribiere (RCG 标准变体)
    options.linesearch = @linesearch_adaptive; % 自适应线搜索 (通常更快)
%     checkgradient(problem);
    % 5. 启动 RCG 求解
    % fprintf('  > R-Step (RCG)...');
    [R_opt, cost_min, info] = conjugategradient(problem, R_init, options);
    % fprintf(' Done. Cost: %.4f\n', cost_min);
    
    % 6. 计算优化后的指标 (方便外层查看)
    [~, sinr_vals] = calc_objective(R_opt, H_func_data, W, sim_param);
    min_sinr_val = min(sinr_vals);
end