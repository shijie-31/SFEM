function showmesh(node,elem,options)
% SHOWMESH 显示二维/三维网格并应用红-白-蓝渐变色方案
% 每个单元随机分配渐变色带中的颜色，保持科学可视化规范

% 参数处理
if nargin==2, options.FaceAlpha = 0.4; end
if nargin==3 && ~isfield(options,'FaceAlpha') 
    options.FaceAlpha = 0.4;
end

% 维度检测
dim = size(node,2);

%% 网格绘制核心逻辑
if ~iscell(elem)
    if size(elem,2)==3  % 三角形网格
        h = patch('Faces', elem, 'Vertices', node);
    elseif size(elem,2)==4 % 四面体网格
        h = tetramesh(elem,node,ones(size(elem,1),1));
    end
else
    if iscell(elem{1}), elem = vertcat(elem{:}); end  % 展平嵌套元胞数组
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % 填充函数
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    h = patch('Faces', tpad, 'Vertices', node);
end

%% 三维可视化设置
if dim==3
    view(3); 
    set(h,'FaceAlpha',options.FaceAlpha); % 透明度设置
end

%% 红-白-蓝渐变色生成 (科学优化版)
numColors = 1024;  % 颜色分辨率
mid = floor(numColors/2);

% 基准色定义 (符合Nature期刊配色规范)
red = [0.8, 0.1, 0.1];   % 增强红色可区分性
white = [1.0, 1.0, 1.0]; % 纯白中点
blue = [0.1, 0.1, 0.9];  % 优化蓝色通道

% 红->白渐变段 (加入gamma校正)
gamma_red = 0.7;
t_red = linspace(0,1,mid)'.^gamma_red;
red_white = [red(1) + t_red*(white(1)-red(1)),...
            red(2) + t_red*(white(2)-red(2)),...
            red(3) + t_red*(white(3)-red(3))];

% 白->蓝渐变段 (sigmoid过渡)
t_blue = 1./(1+exp(-6*(linspace(-0.5,0.5,numColors-mid)')));
white_blue = [white(1) + t_blue*(blue(1)-white(1)),...
             white(2) + t_blue*(blue(2)-white(2)),...
             white(3) + t_blue*(blue(3)-white(3))];

% 合并为完整colormap
coolwarm_reversed = [red_white; white_blue];

%% 颜色分配系统
num_elems = size(elem, 1);  % 单元总数

% 生成非均匀随机分布（避免颜色聚集）
rand_samples = rand(num_elems,1);
color_indices = ceil(numColors * (rand_samples.^2)); % 平方分布偏向深色区域

% 单元颜色矩阵
element_colors = coolwarm_reversed(color_indices, :);

%% 顶点颜色映射
if iscell(elem)
    num_vertices = size(node, 1);
    vertex_colors = zeros(num_vertices, 3);
    for i = 1:num_elems
        vertex_indices = elem{i};
        vertex_colors(vertex_indices, :) = repmat(element_colors(i,:),...
            length(vertex_indices), 1);
    end
else
    num_vertices = size(node, 1);
    vertex_colors = zeros(num_vertices, 3);
    for i = 1:num_elems
        vertex_indices = elem(i,:);
        vertex_colors(vertex_indices, :) = repmat(element_colors(i,:),...
            length(vertex_indices), 1);
    end
end

%% 图形属性设置
set(h,...
    'FaceVertexCData', vertex_colors,... % 顶点颜色数据
    'FaceColor', 'flat',...              % 平坦着色模式
    'EdgeColor', 'k',...                 % 黑色网格线
    'LineWidth', 0.5);                   % 线宽调整

% 默认面颜色设置
if isfield(options,'facecolor') 
    facecolor = options.facecolor;
else
    facecolor = [0.9 0.9 0.9]; % 浅灰背景
end
set(gcf,'Color',facecolor);

%% 视图优化
axis equal tight;
sh = 0.02; % 边界缓冲系数
if dim==2
    xlim([min(node(:,1))-sh, max(node(:,1))+sh])
    ylim([min(node(:,2))-sh, max(node(:,2))+sh])
else
    zlim([min(node(:,3))-sh, max(node(:,3))+sh])
end
axis off; % 关闭坐标轴

%% 色标添加（可选）
% colorbar('Ticks',[0 1],...
%          'TickLabels',{'Hot (Red)','Cold (Blue)'},...
%          'FontSize',11);
end
