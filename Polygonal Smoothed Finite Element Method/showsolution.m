function showsolution(node,elem,u,varargin)
%Showsolution displays the solution corresponding to a mesh given by [node,elem] in 2-D.
%
% Copyright (C) Terence Yu.

data = node;
if ~iscell(elem)
    patch('Faces', elem,...
        'Vertices', data,...
        'FaceColor', 'interp',...
        'CData', u);
else
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    patch('Faces', tpad,...
        'Vertices', data,...
        'FaceColor', 'interp',...
        'CData', u);
end
axis equal; 
sh = 0.0;
xlim([min(node(:,1)) - sh, max(node(:,1)) + sh])
ylim([min(node(:,2)) - sh, max(node(:,2)) + sh])

xlabel('x'); ylabel('y'); 
axis equal; 
sh = 0.0;
xlim([min(node(:,1)) - sh, max(node(:,1)) + sh])
ylim([min(node(:,2)) - sh, max(node(:,2)) + sh])

xlabel('x'); ylabel('y'); 
%colormap('parula');  % 适合科学绘图的内置色彩图
numColors = 1024;  % 色彩图的分辨率
mid = floor(numColors/2);

% 从蓝色到白色
% 定义新的基准色（红-白-蓝）
red = [0.7, 0, 0];      % 深红起始
white = [1.0, 1.0, 1.0];% 纯白中点
blue = [0, 0, 0.7];     % 深蓝结束

blue_to_white = [linspace(blue(1), white(1), mid)', ...
                linspace(blue(2), white(2), mid)', ...
                linspace(blue(3), white(3), mid)'];

% 从白色到红色
white_to_red = [linspace(white(1), red(1), numColors - mid)', ...
               linspace(white(2), red(2), numColors - mid)', ...
               linspace(white(3), red(3), numColors - mid)'];

% 合并两个部分
coolwarm = [blue_to_white; white_to_red];

colormap(coolwarm);
%camlight('headlight');
%lighting gouraud;  % 使用 Gouraud 光照实现平滑阴影



