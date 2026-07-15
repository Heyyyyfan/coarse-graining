function out = biofilm_fractal(stlPath)
% biofilm_fractal 估计生物膜 STL 表面的分形维数，并给出网格类型分析
%
% 用法:
%   out = biofilm_fractal('3xfood-1- 5 days_c5.stl');
%
% 输出 out 结构体包含：
%   .isWatertight   是否水密
%   .nV, .nF        顶点/三角形数量
%   .bbox           归一化前的包围盒 [min; max]
%   .scaleInfo      归一化缩放和平移
%   .boxcount       盒计数的尺度、计数与线性回归结果
%   .D_box          盒计数分形维数估计
%   .coarse         粗粒度尺度/面积等（同时含归一化与绝对值列）
%   .D_coarse       基于粗粒度的维数估计（工程近似）
%
% 备注：
%  - STL 若非水密，体素填充体积不可靠，但盒计数对“表面相交体素”的统计仍可用。
%  - 粗粒度模块以体素形态学开闭代替论文中的“嵌套表面融合”，是工程近似实现。
%  - 本版仅“在输出阶段”把归一化值还原为原模型坐标下的绝对值，不改变任何计算流程。

%% 1) 读取 STL
[V,F] = localReadSTL(stlPath);           % V: Nx3, F: Mx3
nV = size(V,1); nF = size(F,1);

% 记录原始包围盒
bbMin = min(V,[],1); bbMax = max(V,[],1);
bbox = [bbMin; bbMax];

% 2) 网格归一化到单位立方体 [0,1]^3，便于跨尺度
scale = max(bbMax - bbMin);              % ← 原模型“特征长度”= 最大边长
Vn = (V - bbMin) ./ scale;               % 平移+等比缩放到 [0,1]^3
shift = bbMin; s = scale;

% 3) 网格水密性检测（每条无向边是否恰好被两片面共享）
isWatertight = localWatertight(F);

% 4) 盒计数（对表面进行体素覆盖统计）
resList = round(2.^(4:9));               % 16,32,64,128,256,512（可按显存/时间调）
[epsList, Nboxes] = localBoxCountSurface(Vn, F, resList);

% 线性回归：log N(ε) = D * log(1/ε) + c
x = log(1./epsList(:));
y = log(Nboxes(:) + eps(1));             % 避免 log(0)
P = polyfit(x,y,1);
D_box = P(1);

%% 5) 粗粒度：体素化-> 形态学“熔化”-> 多尺度 At 与 Ae
doCoarse = true;
coarse = struct([]);
D_coarse = NaN;
baseRes = 512;                           % 显存不够可降到 384
[shell0, vx] = localVoxelizeShell(Vn, F, baseRes);
shell0 = imdilate(shell0, strel('cube',2));  % 轻度加厚，避免断裂

% —— 用“边界背景标记法”构造实心体 ——
bg = ~shell0;
CCb = bwconncomp(bg, 26);
[nX,nY,nZ] = size(bg);
isBorder = false(CCb.NumObjects,1);
for k = 1:CCb.NumObjects
    idx = CCb.PixelIdxList{k};
    [x1,y1,z1] = ind2sub([nX,nY,nZ], idx);
    if any(x1==1 | x1==nX | y1==1 | y1==nY | z1==1 | z1==nZ)
        isBorder(k) = true;
    end
end
bg_border = false(size(bg));
for k = find(isBorder).'
    bg_border(CCb.PixelIdxList{k}) = true;
end
bg_interior = bg & ~bg_border;           % 内部空洞背景
solid0 = shell0 | bg_interior;           % “壳” ∪ “内部空洞” = 实心体
solid0 = imdilate(solid0, strel('cube',1));  % 再轻度加厚

% —— 方案A：在粗粒度之前，算一份“voxel 版”的原始表面积 A0_vox（归一化） ——
pad0 = 2;
solidP0 = padarray(solid0, [pad0 pad0 pad0], 0, 'both');
iso0 = isosurface(solidP0, 0.5);
if isempty(iso0.vertices)
    A0_vox = NaN;
else
    Viso0 = (iso0.vertices - pad0) * vx;
    A0_vox = sum(localTriArea(Viso0, iso0.faces));   % 归一化坐标下的体素等值面面积
end
% —— 方案A结束 ——

% —— 粗粒度尺度（归一化长度）——
seVox  = [4 6 8 12 16 24 32 48]';        % 起点从 4，点不够再加 48
lambda = (2*seVox + 1) * vx;             % 归一化尺度

At = nan(numel(seVox),1);                % 外表面积（归一化）
Ae = nan(numel(seVox),1);                % 包络面积（归一化）
Vg = nan(numel(seVox),1);                % “灰质体积”（归一化）
Tth = nan(numel(seVox),1);               % 厚度（归一化）= Vg/At

for ii = 1:numel(seVox)
    r  = seVox(ii);
    se = strel('cube', max(1,2*r+1));

    % 在“实心体”上做粗粒化：先闭再开（先填小凹，再去小毛刺），都很温和
    solid = imclose(solid0, se);
    solid = imopen (solid,  strel('cube',3));

    % 占用率守卫：跳过接近全 0 / 全 1 的尺度
    occ = nnz(solid)/numel(solid);
    if occ < 0.005 || occ > 0.995
        At(ii)=NaN; Ae(ii)=NaN; Vg(ii)=NaN; Tth(ii)=NaN; continue;
    end

    % 外侧补零，保证存在外部背景
    pad = 2;
    solidP = padarray(solid, [pad pad pad], 0, 'both');

    % 只取“外表面”
    iso = isosurface(solidP, 0.5);
    if isempty(iso.vertices)
        At(ii)=NaN; Ae(ii)=NaN; Vg(ii)=NaN; Tth(ii)=NaN; continue;
    end

    Viso = (iso.vertices - pad) * vx;
    Atri = localTriArea(Viso, iso.faces);    At(ii) = sum(Atri);

    % —— Ae(λ)：对外部的内部距离剥壳（自适应稳健版，归一化值）——
    bgP = ~solidP;
    bgInterior = imclearborder(bgP, 26);     % 内部空洞
    bgExternal = bgP & ~bgInterior;          % 外部背景

    D_to_ext = bwdist(bgExternal);           % 欧氏距离（体素）
    rvox_base = max(1, round(0.1*r));        % 比 r 小，避免过度
    frac_min = 0.03;                         % 至少保留 3% 体素
    nSolid = nnz(solidP);
    rvox_try = rvox_base;

    core = D_to_ext > rvox_try;
    tries = 0;
    while nnz(core) < frac_min*nSolid && rvox_try > 1 && tries < 3
        rvox_try = max(1, rvox_try - 1);
        core = D_to_ext > rvox_try;
        tries = tries + 1;
    end

    if nnz(core) < max(100, 0.005*nSolid)
        Ae(ii) = NaN;
    else
        isoE = isosurface(core, 0.5);
        if isempty(isoE.vertices)
            Ae(ii) = NaN;
        else
            Venv = (isoE.vertices - pad) * vx;
            Ae_val = sum(localTriArea(Venv, isoE.faces));
            % 数值护栏：Ae ≤ At；并强制随尺度不增
            Ae_val = min(Ae_val, At(ii));
            if ii > 1 && isfinite(Ae(ii-1)), Ae_val = min(Ae_val, Ae(ii-1)); end
            Ae(ii) = Ae_val;
        end
    end

    % —— 体积与厚度（归一化值）——
    voxVol = vx^3;
    Vg(ii) = nnz(solid) * voxVol;            % “灰质体积”工程近似
    if isfinite(At(ii)) && At(ii) > 0
        Tth(ii) = Vg(ii) / At(ii);           % 厚度（归一化）
    else
        Tth(ii) = NaN;
    end
end

% 自检输出（归一化值）
disp(table(seVox(:), lambda(:), At(:), Ae(:), ...
    'VariableNames', {'seVox','lambda_norm','At_norm','Ae_norm'}));
fprintf('有效点数(At>0) = %d\n', nnz(isfinite(At)&At>0));

% === 仅改变“输出”：把归一化值还原为原模型坐标下的绝对值 ===
lenScale  = s;             % 长度缩放因子（= 原模型最大边长）
areaScale = lenScale^2;    % 面积还原
volScale  = lenScale^3;    % 体积还原

lambda_abs = lambda(:) * lenScale;
At_abs     = At(:)     * areaScale;
Ae_abs     = Ae(:)     * areaScale;
Vg_abs     = Vg(:)     * volScale;
T_abs      = Tth(:)    * lenScale;     % T = V/A → × s

% 方案A：voxel 版原始面积的绝对值 / 归一化值
A0_vox_norm = A0_vox;                 % 归一化空间里的面积
A0_vox_abs  = A0_vox * areaScale;     % 换回原始坐标的面积

% —— 存表（同时含归一化与绝对值，便于论文引用）——
coarse = table( ...
    lambda(:), lambda_abs(:), ...
    At(:),     At_abs(:), ...
    Ae(:),     Ae_abs(:), ...
    Vg(:),     Vg_abs(:), ...
    Tth(:),    T_abs(:), ...
    'VariableNames', { ...
        'lambda_norm','lambda_abs', ...
        'At_norm','At_abs', ...
        'Ae_norm','Ae_abs', ...
        'Vg_norm','Vg_abs', ...
        'T_norm','T_abs' ...
    });
writetable(coarse, 'coarse_lambda_At_Ae_T.csv');

%% 6) 打印摘要
fprintf('== STL 网格信息 ==\n');
fprintf(' 顶点: %d, 三角形: %d\n', nV, nF);
fprintf(' 包围盒原始尺寸: [%.3f %.3f %.3f]\n', (bbMax-bbMin));
fprintf(' 网格是否水密: %d (1=是,0=否)\n', isWatertight);
fprintf(' 盒计数分形维数 D_box = %.4f\n', D_box);

%% 7) 输出结构体（包含两套值）
out = struct();
out.isWatertight = isWatertight;
out.nV = nV; out.nF = nF;
out.bbox = bbox;
out.scaleInfo.shift = shift; out.scaleInfo.scale = s;
out.boxcount.resList = resList;
out.boxcount.epsList = epsList;
out.boxcount.Nboxes = Nboxes;
out.D_box = D_box;

% 方案A：附加一份“体素版原始表面积”
out.A0_vox_norm = A0_vox_norm;   % 归一化坐标下的 voxel 表面积
out.A0_vox_abs  = A0_vox_abs;    % 原始物理坐标下的 voxel 表面积

if doCoarse
    out.coarse = coarse;                 % 同时有 *_norm 与 *_abs
end

%% 8) 分形维（仍用归一化 At 拟合，不受还原影响）
% 仅保留有效点，并按尺度从小到大排序（归一化）
ok = isfinite(lambda) & isfinite(At) & At>0;
lam_fit = lambda(ok);
A_fit   = At(ok);
[lam_fit,ord] = sort(lam_fit);  A_fit = A_fit(ord);

% 轻度单调不增修正（仅用于稳健拟合；可注释掉）
Amono = A_fit;
for i = 2:numel(Amono)
    if Amono(i) > Amono(i-1), Amono(i) = Amono(i-1); end
end

% 回归 D = 2 + slope(log At, log(1/λ))（用归一化数值，斜率与是否还原无关）
xfit = log(1./lam_fit);
yfit = log(Amono);
if numel(xfit) >= 3
    p = polyfit(xfit, yfit, 1);
    m = max(0, min(1, p(1)));           % 斜率守卫到 [0,1]，对应 D ∈ [2,3]
    D_coarse = 2 + m;
else
    D_coarse = NaN;
end

fprintf(' 粗粒度(工程近似)分形维数 D_coarse = %.4f\n', D_coarse);
out.D_coarse = D_coarse;

%% 9) 可视化（仍画归一化尺度 vs 归一化面积；如需看绝对值可改 y 数据为 At_abs）
try
    figure('Color','w');
    plot(lam_fit, A_fit, 'o-', 'LineWidth', 2, 'MarkerSize', 6); grid on;
    xlabel('\lambda (normalized)'); ylabel('A_t(\lambda) (normalized)');
    title(sprintf('Coarse-graining: A_t vs. \\lambda   (D_{coarse}=%.4f)', D_coarse));
    exportgraphics(gcf, 'coarse_At_vs_scale.png', 'Resolution', 300);
catch
end

% 另存回归点（归一化）
T_points = table(lam_fit, A_fit, log(1./lam_fit), log(A_fit), ...
    'VariableNames', {'lambda_norm','At_norm','log_inv_lambda','log_At_norm'});
writetable(T_points, 'coarse_At_vs_scale_points.csv');

end

%% ========= 辅助函数 =========
function [V,F] = localReadSTL(p)
% 兼容多种 stlread 返回类型：triangulation / struct / [F,V]
    triedAscii = false;
    try
        tr = stlread(p);
        if isa(tr, 'triangulation')
            V = tr.Points;
            F = tr.ConnectivityList;
        elseif isstruct(tr)
            if isfield(tr,'Vertices') && isfield(tr,'Faces')
                V = tr.Vertices;  F = tr.Faces;
            elseif isfield(tr,'vertices') && isfield(tr,'faces')
                V = tr.vertices;  F = tr.faces;
            else
                error('未知的 struct 形式。');
            end
        elseif isnumeric(tr)
            try
                [F,V] = stlread(p); % 再试双输出
            catch
                error('stlread 返回的数值类型无法解析。');
            end
        else
            error('不支持的 stlread 返回类型。');
        end
    catch
        % 退回到简易 ASCII STL 解析（若是二进制 STL 会失败）
        triedAscii = true;
        fid = fopen(p,'r');
        if fid<0, error('无法打开 STL 文件。'); end
        C = textscan(fid,'%s','Delimiter','\n'); fclose(fid);
        C = string(C{1});
        i = contains(lower(C),'vertex');
        if ~any(i)
            error(['无法解析 STL。若为二进制 STL，请使用 MATLAB 自带 stlread，' ...
                   '或在 Meshlab/Blender 中另存为 ASCII STL 再试。']);
        end
        Vlist = sscanf(strjoin(extractAfter(C(i),"vertex "), newline), ...
                       '%f %f %f', [3, Inf])';
        if mod(size(Vlist,1),3)~=0
            error('ASCII STL 顶点数不是 3 的倍数，文件可能损坏或为二进制。');
        end
        V = Vlist;
        F = reshape(1:size(Vlist,1), 3, [])';
    end

    % 统一为 double 索引
    V = double(V);
    F = double(F);

    % 去重顶点并重映射
    [V,~,ix] = unique(V, 'rows', 'stable');
    F = reshape(ix(F), size(F));

    % 基本健壮性检查
    if size(F,2) ~= 3
        % 有些实现会返回非三角网格（极少见），尝试三角化
        try
            TR = triangulation(F, V);
            F = TR.ConnectivityList;
        catch
            error('网格不是三角形面片，且无法自动三角化。');
        end
    end

    if isempty(V) || isempty(F)
        if triedAscii
            error('未能从 ASCII STL 解析出顶点/面。');
        else
            error('stlread 解析失败：请确认文件是否有效 STL。');
        end
    end
end

function tf = localWatertight(F)
% 判断三角网格是否水密：每条无向边是否恰好被两片面共享
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = sort(E,2);
    [u,~,ic] = unique(E,'rows','stable');
    ic = double(ic(:));
    nEdges = size(u,1);
    cnt = accumarray(ic, 1, [nEdges 1], @sum);
    tf = all(cnt == 2);
end

function [epsList, Nboxes] = localBoxCountSurface(V,F,resList)
% 更精确的表面盒计数：以体素中心为采样点，判断与任意三角形的最近距离是否 < 半个对角线
    epsList = 1./resList;
    Nboxes  = zeros(size(resList));

    % 预计算每个三角形的法向和边向量
    T1 = V(F(:,2),:) - V(F(:,1),:);
    T2 = V(F(:,3),:) - V(F(:,1),:);

    for k = 1:numel(resList)
        R = resList(k);
        eps = 1/R;                            % 盒子边长（归一化）
        rchk = sqrt(3)*eps/2;                 % 体素中心到表面判定阈值（体素外接球半径）

        % 体素中心坐标网格
        c = ( (0:R-1) + 0.5 ) * eps;          % 0.5 偏置是体素中心
        [X,Y,Z] = ndgrid(c,c,c);
        occ = false(R,R,R);

        % 粗筛：为每个三角形算 AABB 对应的体素范围
        Vi = max(1, min(R, floor(V.*R) + 1));
        %#ok<NASGU>
        for t = 1:size(F,1)
            tri = V(F(t,:),:);
            mn = max([1 1 1], floor(min(tri,[],1).*R) );
            mx = min([R R R], ceil (max(tri,[],1).*R) );
            if any(mx<mn), continue; end

            % 拉出这个小块的体素中心
            x = X(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
            y = Y(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
            z = Z(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
            P = [x(:), y(:), z(:)];                           % Q 个点

            % 近似：点到三角面所在平面的距离 + 投影重心坐标判定
            p0 = tri(1,:); e0 = T1(t,:); e1 = T2(t,:);
            n = cross(e0, e1); n = n / norm(n);
            distPlane = abs((P - p0) * n');                   % Qx1

            Pproj = P - distPlane.*n;
            u = Pproj - p0;
            M = [e0; e1]'; st = (M\u')'; s = st(:,1); t2 = st(:,2);
            inside = s>=0 & t2>=0 & s+t2<=1;
            hit = inside & (distPlane <= rchk);

            blk = false(size(x));
            blk(hit) = true;
            occ(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3)) = ...
                occ(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3)) | blk;
        end
        Nboxes(k) = nnz(occ);
    end
end

function [vol, vx] = localVoxelizeShell(V, F, R)
% 基于三角形 AABB 的“壳层体素化”，得到体素逻辑阵（占用表示靠近表面）
% 归一化空间下体素尺寸 vx。
    vol = false(R,R,R);
    Vi = max(1, min(R, floor(V.*R) + 1));
    pad = 1;  % 给壳层一点厚度（±1 体素）
    for t = 1:size(F,1)
        tri = Vi(F(t,:),:);
        mn = max([1 1 1], min(tri,[],1)-pad);
        mx = min([R R R], max(tri,[],1)+pad);
        vol(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3)) = true;
    end
    vx = 1/R; % 归一化下的体素边长
end

function A = localTriArea(V,F)
% 计算三角面片面积
    v1 = V(F(:,2),:) - V(F(:,1),:);
    v2 = V(F(:,3),:) - V(F(:,1),:);
    A = 0.5*vecnorm(cross(v1,v2,2),2,2);
end
