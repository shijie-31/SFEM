function phi_corrected = mvc(v, x)
% v: 多边形顶点坐标 (n x 2)
% x: 目标点坐标 (1 x 2)
n = size(v, 1); % 顶点数
phi = zeros(n, 1); % 初始化形函数
w = zeros(n, 1); % 初始化权重
gamma = 1; % 调参数，适应凹顶点
% 计算每条边的角度和权重
for i = 1:n
% 获取相邻顶点的索引
im1 = mod(i-2, n) + 1;
ip1 = mod(i, n) + 1;
% 计算相邻顶点到目标点的向量
vi = v(i, :) - x;
vim1 = v(im1, :) - x;
vip1 = v(ip1, :) - x;
% 规范化向量
vi = vi / norm(vi);
vim1 = vim1 / norm(vim1);
vip1 = vip1 / norm(vip1);
% 计算相邻向量之间的夹角
theta_im1 = acos(max(-1, min(1, dot(vim1, vi))));
theta_ip1 = acos(max(-1, min(1, dot(vip1, vi))));
% 计算当前顶点的权重
w(i) = (tan(theta_im1 / 2) + tan(theta_ip1 / 2)) / (norm(v(i, :) - x)^gamma);
end
% 归一化权重，得到形函数
phi_corrected = w / sum(w);
end