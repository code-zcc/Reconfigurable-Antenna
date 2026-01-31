function R_out = rotate_frame_z_to_r3(r3_in)
%   输入: 
%       r3_in: [3 x 1] 向量 或 [3 x M] 矩阵 (每一列是一个 r3)
%   输出:
%       R_out: [3 x 3] 矩阵 或 [3 x 3 x M] 多维数组
%
%   构造规则:
%       使用最小旋转 (Minimal Rotation) 将 Z轴 [0;0;1] 映射到 r3。
%       生成的 R 满足: R(:,3) == r3 (归一化后)。
%       r1, r2 由右手定则自然生成，无额外 Roll。

    [~, M] = size(r3_in);
    R_out = zeros(3, 3, M);
    z_ref = [0; 0; 1];

    for m = 1:M
        target_z = r3_in(:, m);
        

        target_z = target_z / norm(target_z);
        

        v = cross(z_ref, target_z);
        s = norm(v);       % sin(theta)
        c = dot(z_ref, target_z); % cos(theta)

        if s < 1e-6
            if c > 0

                R_out(:, :, m) = eye(3);
            else

                R_out(:, :, m) = diag([1, -1, -1]);
            end
        else

            vx = [0,      -v(3),  v(2); ...
                  v(3),   0,     -v(1); ...
                 -v(2),   v(1),   0   ];
             
            R_out(:, :, m) = eye(3) + vx + vx^2 * ((1 - c) / s^2);
        end
    end
    

    if M == 1
        R_out = squeeze(R_out);
    end
end