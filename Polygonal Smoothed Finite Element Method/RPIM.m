function [phi, dphi] = rpim_shape_function(v, x)
% 计算 RPIM（径向点插值法）下的形函数值和梯度
% 输入参数：
% v : [x1 y1; x2 y2; x3 y3]，三角形单元的三个顶点坐标
% x : [x, y]，待计算形函数的点的坐标
% 输出参数：
% phi : 形函数值向量 [phi1; phi2; phi3]
% dphi : 形函数的梯度矩阵 [dphi1_dx, dphi1_dy; dphi2_dx, dphi2_dy; dphi3_dx, dphi3_dy]

% 参数设置
ac = 0.35;
q = 1.03;

% 计算 dc，节点间最大距离
dc1 = sqrt((v(1,1) - v(2,1))^2 + (v(1,2) - v(2,2))^2);
dc2 = sqrt((v(2,1) - v(3,1))^2 + (v(2,2) - v(3,2))^2);
dc3 = sqrt((v(3,1) - v(1,1))^2 + (v(3,2) - v(1,2))^2);
dc = max([dc1, dc2, dc3]);

% 定义径向基函数及其导数
rbf = @(xx, xi, yy, yi) ((xx - xi)^2 + (yy - yi)^2 + (ac * dc)^2)^q;
drbf = @(xx, xi, yy, yi) 2 * q * ((xx - xi)^2 + (yy - yi)^2 + (ac * dc)^2)^(q - 1);

% 节点数
wl = size(v, 1);

% 构建 RBF 矩阵 Rq 和多项式矩阵 Pm
Rq = zeros(wl, wl);
Pm = zeros(wl, 3);
for k = 1:wl
    for j = 1:wl
        Rq(k, j) = rbf(v(k,1), v(j,1), v(k,2), v(j,2));
    end
    Pm(k, :) = [1, v(k,1), v(k,2)];
end

% 计算形函数系数 Sa 和 Sb
invRq = inv(Rq);
Sb = inv(Pm' * invRq * Pm) * Pm' * invRq;
Sa = invRq - invRq * Pm * Sb;

% 在点 x 处计算形函数值
Rqsg = zeros(1, wl);
for w = 1:wl
    Rqsg(w) = rbf(x(1), v(w,1), x(2), v(w,2));
end
Pmsg = [1, x(1), x(2)];
phi = (Rqsg * Sa + Pmsg * Sb)';

% 计算形函数的梯度
dRqsg = zeros(wl, 2);
for w = 1:wl
    dx = x(1) - v(w,1);
    dy = x(2) - v(w,2);
    common_term = drbf(x(1), v(w,1), x(2), v(w,2));
    dRqsg(w,1) = common_term * dx;
    dRqsg(w,2) = common_term * dy;
end

dPmsg = [0, 1, 0; 0, 0, 1];
dphi = zeros(wl, 2);
for dim = 1:2
    dphi(:, dim) = (dRqsg(:, dim)' * Sa + dPmsg(dim, :) * Sb)';
end

end