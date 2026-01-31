function R = get_rotation_matrix_3d(roll, pitch, yaw)
% 功能:
%   根据 Z-Y-X (Yaw-Pitch-Roll) 旋转顺序生成旋转矩阵 R。
% 输入:
%   roll  - [Scalar] 滚转角 (绕 X 轴), rad
%   pitch - [Scalar] 俯仰角 (绕 Y 轴), rad
%   yaw   - [Scalar] 偏航角 (绕 Z 轴), rad
% 输出:
%   R     - [3 x 3] 旋转矩阵
    Rx = [1 0 0; 
          0 cos(roll) -sin(roll); 
          0 sin(roll) cos(roll)];
      

    Ry = [cos(pitch) 0 sin(pitch); 
          0 1 0; 
          -sin(pitch) 0 cos(pitch)];
      

    Rz = [cos(yaw) -sin(yaw) 0; 
          sin(yaw) cos(yaw) 0; 
          0 0 1];
      
    R = Rz * Ry * Rx;
end