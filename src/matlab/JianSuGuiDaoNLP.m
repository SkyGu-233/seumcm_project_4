function rocket_optimization_solution
    % 参数设置
    m0 = 1000;     % 初始质量 (kg)
    v_e = 3000;    % 排气速度 (m/s)
    g = 9.8;       % 重力加速度 (m/s^2)
    v0 = 0;        % 初始y方向速度 (m/s)
    x_desired = 1000; % 目标位移 (m)
    T_min = 1500;  % 最小推力 (N)
    T_max = 7500;  % 最大推力 (N)
    
    % 初始猜测 [T_x, T_y, t]
    % T_x0 = 4000, T_y0 = 0 (假设推力主要在x方向), t0=100秒
    X0 = [4000; 0; 100];
    
    % 变量边界
    lb = [-T_max, -T_max, 0];   % T_x, T_y, t 的下界
    ub = [T_max, T_max, Inf];   % T_x, T_y, t 的上界
    
    % 优化选项
    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'interior-point', ... % 内点法处理非线性约束
        'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 1000, ...
        'StepTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6);
    
    % 调用fmincon求解
    [X_opt, fval, exitflag] = fmincon(@(X)objective(X), X0, [], [], [], [], lb, ub, ...
        @(X)constraints(X, m0, v_e, g, v0, x_desired, T_min, T_max), options);
    
    % 输出结果
    if exitflag > 0
        fprintf('优化成功!\n');
        fprintf('最优解:\n');
        fprintf('T_x = %.4f N\n', X_opt(1));
        fprintf('T_y = %.4f N\n', X_opt(2));
        fprintf('t = %.4f s\n', X_opt(3));
        
        % 计算最优目标函数值 (最大化 m1)
        T_opt = sqrt(X_opt(1)^2 + X_opt(2)^2);
        v_m_opt = T_opt / v_e;
        m1_opt = m0 - v_m_opt * X_opt(3);
        fprintf('最优 m1 = %.4f kg\n', m1_opt);
    else
        fprintf('优化未收敛到解。\n');
    end
end

% 目标函数 (最小化 v_m * t)
function f = objective(X)
    T_x = X(1);
    T_y = X(2);
    t = X(3);
    T = sqrt(T_x^2 + T_y^2);
    v_m = T / v_e;  % v_e 在外部参数中定义
    f = v_m * t;    % 等价于最小化 v_m * t
end

% 非线性约束
function [c, ceq] = constraints(X, m0, v_e, g, v0, x_desired, T_min, T_max)
    T_x = X(1);
    T_y = X(2);
    t = X(3);
    
    % 计算推力大小
    T = sqrt(T_x^2 + T_y^2);
    v_m = T / v_e;
    
    % 1. 推力约束 (不等式约束)
    c1 = T_min - T;   % T >= T_min
    c2 = T - T_max;   % T <= T_max
    c3 = v_m * t - m0; % m0 - v_m*t >= 0 (质量非负)
    c = [c1; c2; c3];
    
    % 2. 位移和速度方程 (等式约束)
    % 位移 Δx 方程
    term1 = (T_x / v_m) * log(m0) * t;
    term2 = 0.5 * g * t^2;
    term3 = (T_x / v_m) * ((t - m0 / v_m) * log(m0 - v_m * t) - t);
    Delta_x = term1 + term2 - term3;
    
    % 速度方程
    x_dot = (T_x / v_m) * log(m0) + g * t - (T_x / v_m) * log(m0 - v_m * t);
    y_dot = v0 + (T_y / v_m) * log(m0) - (T_y / v_m) * log(m0 - v_m * t);
    v1 = sqrt(x_dot^2 + y_dot^2);
    
    % 等式约束: Δx = x_desired, v1 无约束 (仅占位)
    ceq = [Delta_x - x_desired; 0]; % 第二个等式约束设为0 (无实际约束)
end