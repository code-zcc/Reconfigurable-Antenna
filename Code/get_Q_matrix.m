function Q = get_Q_matrix(R_rx, e_theta, e_phi)
% 功能:
%   计算从波传播坐标系(e_theta, e_phi)到接收天线本地极化端口的投影矩阵。
% 输入:
%   R_rx    - [3 x 3] 接收天线的旋转矩阵 (定义了天线姿态)
%   e_theta - [3 x 1] 到达波的 theta 分量基向量
%   e_phi   - [3 x 1] 到达波的 phi 分量基向量
% 输出:
%   Q       - [2 x 2] 投影矩阵
    u_loc_1 = [1; 0; 0];
    u_loc_2 = [0; 1; 0];
    

    u_glob_1 = R_rx * u_loc_1;
    u_glob_2 = R_rx * u_loc_2;
    

    Q = [dot(u_glob_1, e_theta), dot(u_glob_1, e_phi);
         dot(u_glob_2, e_theta), dot(u_glob_2, e_phi)];
end