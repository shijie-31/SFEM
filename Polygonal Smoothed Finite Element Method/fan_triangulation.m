function [V,T,tri_count] = fan_triangulation(e_node, nod)
% FAN_TRIANGULATION  质心‑顶点扇形三角剖分
%
% [V,T] = fan_triangulation(e_node)
%
% INPUT
%   e_node : n×2 多边形顶点坐标（按顺/逆时针）
%
% OUTPUT
%   V : (n+1)×2 所有节点坐标，最后一行是多边形质心
%   T : n×3    三角形索引，每行为 [centerIdx, vi, v(i+1)]
% ------------------------------------------------------------
% 最后更新: 2025‑05‑08

% ---------- 1. “坐标平均”求中心 ----------
center = mean(e_node, 1);      % 简单平均而非面积加权

% ---------- 2. 构造节点与三角形 ----------
V = [e_node; center];          % 节点表
n = size(e_node, 1);
centerIdx = n + 1;
T = [repmat(centerIdx, n, 1), (1:n)', [2:n, 1]'];  % 扇形三角单元

tri_count = size(T, 1);
end
