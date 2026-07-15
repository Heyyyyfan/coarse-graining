import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import tifffile
from stl import mesh # 用于生成STL文件

plt.rcParams.update({'font.size': 18})

def process_biofilm_image(image_path):
    """
    处理细菌生物膜图像，计算厚度并三维可视化

    参数:
        image_path: 输入的TIFF图像路径

    返回:
        thickness_map: 厚度矩阵
    """
    # 1. 读取TIFF图像
    image = tifffile.imread(image_path)

    # 确保图像是二维的
    if len(image.shape) > 2:
        image = image[:, :, 0]  # 如果是RGB，取第一个通道

    # 2. 将灰度值转换为OD值 (式4-2)
    # 假设最大灰度值为255（8位图像）
    max_gray = 65530
    image = np.clip(image,1, max_gray)
    od_map = -np.log10((image) / max_gray)  # 加1避免log(0)


    # 3. 将OD值转换为厚度 (式4-3)
    thickness_map = (od_map) * 1.2 - 1.6
    thickness_map[thickness_map > 3] = 0
    thickness_map[thickness_map < 0] = 0
    return thickness_map


def visualize_3d_surface(thickness_map):
    """
    三维可视化厚度图

    参数:
        thickness_map: 厚度矩阵
    """
    # 创建网格
    cols, rows = thickness_map.shape[1], thickness_map.shape[0]
    x = np.linspace(0, X, cols)
    y = np.linspace(0, Y, rows)
    x, y = np.meshgrid(x, y)

    # 创建3D图形
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 绘制表面图
    surf = ax.plot_surface(x, y, thickness_map,
                           cmap=cm.coolwarm,
                           linewidth=0,
                           antialiased=True,
                           rstride=5,
                           cstride=5)

    ax.set_zlim(0, 5)

    # 修改颜色条标签
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Thickness (mm)')

    # 设置标签(单位改为mm)
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Thickness (mm)')

    plt.tight_layout()
    plt.show()


def process_lamda_grid(thickness_map, lamda, X, Y):
    """
    将空间划分为边长为lamda的正方形网格，并遍历每个网格单元

    参数:
        thickness_map: 厚度矩阵
        lamda: 网格边长
        X: X方向总长度
        Y: Y方向总长度

    返回:
        grid_data: 网格数据，包含每个网格的信息
    """


    # 计算网格数量
    x_cells = int(np.ceil(X / lamda / 6))
    y_cells = int(np.ceil(Y / lamda / 6))

    # 获取厚度图尺寸
    rows, cols = thickness_map.shape

    # 计算每个像素代表的实际尺寸
    pixel_size_x = X / cols
    pixel_size_y = Y / rows
    #print(rows, cols)

    # 初始化网格数据存储
    grid_data = np.zeros((y_cells, x_cells))



    # 遍历每个网格单元
    for i in range(y_cells):
        for j in range(x_cells):
            # 计算当前网格的像素范围
            y_start = int(i * lamda * 6 / pixel_size_y)
            y_end = int(min((i + 1) * lamda * 6 / pixel_size_y, rows))

            x_start = int(j * lamda * 6 / pixel_size_x)
            x_end = int(min((j + 1) * lamda * 6 / pixel_size_x, cols))

            z_steps = int(np.ceil(thickness_map[y_start, x_start] / lamda)-1)

            inside_tisu = 0
            if thickness_map[y_start, x_start] > z_steps * lamda:
                inside_tisu  +=  1
            if thickness_map[y_start, x_start] > (z_steps + 1) * lamda:
                inside_tisu  +=  1
            if thickness_map[y_start, x_end - 1] > z_steps * lamda :
                inside_tisu  +=  1
            if thickness_map[y_start, x_end - 1] > (z_steps + 1) * lamda:
                inside_tisu  +=  1
            if thickness_map[y_end - 1, x_start] > z_steps * lamda:
                inside_tisu  +=  1
            if thickness_map[y_end - 1, x_start] > (z_steps + 1) * lamda:
                inside_tisu  +=  1
            if thickness_map[y_end - 1, x_end - 1] > z_steps * lamda:
                inside_tisu  +=  1
            if thickness_map[y_end - 1, x_end - 1] > (z_steps + 1) * lamda :
                inside_tisu  +=  1
            #if (thickness_map[y_start, x_start] < z_steps * lamda and
            #        thickness_map[y_start, x_end - 1] < z_steps * lamda and
            #        thickness_map[y_end - 1, x_start] < z_steps * lamda and
            #        thickness_map[y_end - 1, x_end - 1] < z_steps * lamda):
            #    inside_tisu = 0
            if inside_tisu > 3:
                target_value = (z_steps + 1 ) * lamda
                modified_count = 0  # 计数器
                for y in range(y_start, y_end):
                    for x in range(x_start, x_end):
                        if thickness_map[y, x] < target_value:
                            thickness_map[y, x] = target_value
                            modified_count += 1

                if modified_count > 0:
                    print(f"网格({i},{j}): 目标值={target_value:.3f}")


    return thickness_map
def visualize_2d_heatmap(thickness_map):
    """
    二维可视化厚度热图

    参数:
        thickness_map: 厚度矩阵
    """
    plt.figure(figsize=(10, 8))
    # 显示时设置extent参数使坐标轴显示0-20mm
    plt.imshow(thickness_map, cmap='hot', interpolation='nearest',
               extent=[0, X, 0, Y])
    cbar = plt.colorbar(label='Thickness (mm)')
    plt.xlabel('X (mm)')
    plt.ylabel('Y (mm)')
    plt.title('2D Heatmap of Bacterial Biofilm Thickness')
    plt.show()

def generate_stl(thickness_map, stl_filename, downsample_factor=5):
    """
    生成STL文件
    参数:
        thickness_map: 厚度矩阵(单位：mm)
        stl_filename: 输出的STL文件路径
        downsample_factor: 降采样因子(默认5)
    """
    # 1. 对厚度图进行降采样
    thickness_map = thickness_map[::downsample_factor, ::downsample_factor]

    rows, cols = thickness_map.shape

    # 2. 创建物理坐标网格(0-20mm)
    x_phys = np.linspace(0, X, cols)
    y_phys = np.linspace(0, Y, rows)
    x_phys, y_phys = np.meshgrid(x_phys, y_phys)

    # 3. 创建顶点和面
    vertices = []
    faces = []

    # 4. 使用更高效的网格生成方式
    for i in range(rows - 1):
        for j in range(cols - 1):
            # 当前网格的四个顶点
            v0 = [x_phys[i, j], y_phys[i, j], thickness_map[i, j]]
            v1 = [x_phys[i, j + 1], y_phys[i, j + 1], thickness_map[i, j + 1]]
            v2 = [x_phys[i + 1, j], y_phys[i + 1, j], thickness_map[i + 1, j]]
            v3 = [x_phys[i + 1, j + 1], y_phys[i + 1, j + 1], thickness_map[i + 1, j + 1]]

            # 添加顶点
            current_vertex_count = len(vertices)
            vertices.extend([v0, v1, v2, v3])

            # 添加两个三角形面
            faces.append([current_vertex_count, current_vertex_count + 1, current_vertex_count + 2])
            faces.append([current_vertex_count + 1, current_vertex_count + 2, current_vertex_count + 3])

    # 5. 创建STL网格
    biofilm_mesh = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            biofilm_mesh.vectors[i][j] = vertices[f[j]]

    # 6. 保存STL文件
    biofilm_mesh.save(stl_filename)


# 主程序
if __name__ == "__main__":
    # 输入图像路径
    image_path = "7xfood-4- 5 days_c5.tif"  # 替换为你的图像路径


    X=23.495
    Y=17.6

    lamda = 0.01

    # 处理图像并计算厚度
    thickness_map = process_biofilm_image(image_path)
    thickness_map1 = process_biofilm_image(image_path)
    # 可视化结果
    #visualize_2d_heatmap(thickness_map)
    visualize_3d_surface(thickness_map)


    # 生成STL模型文件
    stl_filename = "7xfood-4- 5 days_c5.stl"
    generate_stl(thickness_map, stl_filename)


    #process_lamda_grid(thickness_map, lamda, X, Y)
    #visualize_2d_heatmap(thickness_map)

    #non_zero_avg = np.mean(thickness_map[thickness_map > 0])
    #print(f"非零区域平均厚度: {non_zero_avg:.3f} mm")
    #non_zero1_avg = np.mean(thickness_map1[thickness_map > 0])
    #print(f"非零区域平均厚度: {non_zero_avg:.3f} mm")

    #visualize_3d_surface(thickness_map)
    # 保存厚度数据
    #np.savetxt("DAY32.csv", thickness_map, delimiter=",")

