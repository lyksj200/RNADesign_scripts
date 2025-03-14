import os
import MDAnalysis as mda
from MDAnalysis.transformations import unwrap, center_in_box, wrap, fit_rot_trans
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysis.analysis import align, pca  # 导入pca模块
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants  # 导入物理常数库
import argparse # 导入 argparse 模块
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA  # 导入 sklearn PCA


# --- 用户需要修改的参数 ---
TEMPERATURE_K = 300  # 模拟温度，单位开尔文 (K)
OUTPUT_FILENAME_PREFIX = 'nucleic_analysis'  # 输出文件名的前缀
FRAME_SAMPLING_INTERVAL = 1 # 轨迹抽帧间隔，每 FRAME_SAMPLING_INTERVAL 帧抽一帧
IMAGE_OUTPUT_DIR_NAME = 'outpngeach'+str(FRAME_SAMPLING_INTERVAL) # 图片输出的文件夹名称
TIME_PER_STEP_NS = 0.02  # 每步模拟的时长，单位纳秒 (ns)

# --- 新增 FEL 分析相关参数 ---
hist_bins_rmsd = 80
hist_bins_rg = 40
hist_bins_pca = 50
gaussian_sigma = 1
figure_size = (24, 6) # 保持原有的 figsize，因为要在一张图中显示四个子图
run_pca_analysis = True # 默认进行 PCA 分析


# --- 参数修改结束 ---

# ---------------------------
#  从参考代码复制的函数定义
# ---------------------------

def load_trajectory(gro_file, xtc_file):
    """加载分子动力学轨迹."""
    return mda.Universe(gro_file, xtc_file)

def calculate_rmsd_rg(universe, selection, ref_universe=None):
    """计算 RMSD 和回转半径 (Rg)."""
    ref_atoms = ref_universe.select_atoms(selection) if ref_universe else universe.select_atoms(selection)
    rmsd_analyzer = RMSD(universe.select_atoms(selection),
                                     ref_atoms,
                                     center=False,
                                     superposition=False).run()
    rmsd_values = rmsd_analyzer.results.rmsd[:, 2]

    rgyr_values = []
    for ts in universe.trajectory:
        rgyr_values.append(universe.select_atoms(selection).radius_of_gyration())
    rgyr_values = np.array(rgyr_values)
    return rmsd_values, rgyr_values

def calculate_pca_sklearn(universe, selection, n_components=2, already_aligned=False):
    """使用 sklearn.decomposition.PCA 执行主成分分析 (PCA)."""
    selected_atoms = universe.select_atoms(selection)
    coords = []
    for ts in universe.trajectory:
        coords.append(selected_atoms.positions.flatten()) # 将每一帧的坐标展平
    data_matrix = np.array(coords)

    pca_sklearn = PCA(n_components=n_components) # 初始化 sklearn PCA
    projections_sklearn = pca_sklearn.fit_transform(data_matrix) # 执行 PCA 拟合和降维
    pc1_sklearn, pc2_sklearn = projections_sklearn[:, 0], projections_sklearn[:, 1] # 获取 PC1 和 PC2 投影

    return pca_sklearn, pc1_sklearn, pc2_sklearn, pca_sklearn.explained_variance_ratio_ # 返回 sklearn PCA 对象, PC1, PC2, 方差解释比例


def calculate_free_energy(x, y, bins_x=80, bins_y=40):
    """计算二维自由能，并允许自定义 x 和 y 轴的 bins 数量，动态设置 bins 范围."""
    # 动态计算 x 和 y 轴的 bins 范围
    x_bins = np.linspace(0.9 * np.min(x), 1.1 * np.max(x), bins_x) # bins_x 个区间
    y_bins = np.linspace(0.9 * np.min(y), 1.1 * np.max(y), bins_y) # bins_y 个区间

    # 使用 numpy.histogram2d 计算二维直方图
    hist, xedge, yedge = np.histogram2d(x, y, bins=[x_bins, y_bins], density=True)
    hist += 1e-12  # 避免零值，防止计算 log 时出现错误
    free_energy = -0.008314 * 300 * np.log(hist)  # 温度300K，单位 KBT
    return free_energy - free_energy.max(), xedge, yedge


def plot_fel(x_edges, y_edges, free_energy, fel_type, variance_explained=None, subplot_ax=None, add_padding=False):
    """
    绘制二维自由能景观图 (FEL 图).

    使用等高线和填充颜色来展示自由能 landscape，颜色从彩虹色到白色表示能量由低到高。

    参数:
    x_edges (numpy.ndarray):  X 轴边缘，由 calculate_free_energy 函数返回
    y_edges (numpy.ndarray):  Y 轴边缘，由 calculate_free_energy 函数返回
    free_energy (numpy.ndarray):  二维自由能矩阵，由 calculate_free_energy 函数返回
    fel_type (str):  自由能图类型 ("PCA" 或 "RMSD-Rg")，用于设置标题和轴标签
    variance_explained (list, 可选):  PCA 的方差解释比例，用于 PCA 图的轴标签，仅当 fel_type="PCA" 时使用，默认为 None
    subplot_ax (matplotlib.axes._axes.Axes, 可选):  matplotlib 子图的 Axes 对象，如果提供，则在指定的子图中绘制，否则在当前 axes 中绘制，默认为 None
    add_padding (bool, 可选): 是否在 FEL 图周围添加留白，默认为 False。 为 PCA 图添加留白通常可以改善视觉效果。
    """
    X, Y = np.meshgrid(x_edges[:-1], y_edges[:-1])
    Z_smooth = gaussian_filter(free_energy.T, sigma=1)  # 对自由能矩阵进行高斯滤波平滑，sigma 值控制平滑程度

    if add_padding: # 判断是否添加留白
        # 设置留白的比例，在图的四周留出空白区域，使等高线和颜色填充不紧贴边缘
        padding_factor = 0.2
        x_padding = (np.max(x_edges) - np.min(x_edges)) * padding_factor
        y_padding = (np.max(y_edges) - np.min(y_edges)) * padding_factor

        # 扩展 x 和 y 轴的范围，用于添加留白
        extended_x_edges = np.linspace(np.min(x_edges) - x_padding, np.max(x_edges) + x_padding, len(x_edges) + int(len(x_edges) * padding_factor * 2))  # 根据 padding_factor 计算扩展的 bins 数量
        extended_y_edges = np.linspace(np.min(y_edges) - y_padding, np.max(y_edges) + y_padding, len(y_edges) + int(len(y_edges) * padding_factor * 2))  # 根据 padding_factor 计算扩展的 bins 数量

        # 创建新的自由能矩阵，并填充最大能量值，用于留白区域
        extended_free_energy = np.full((len(extended_x_edges) - 1, len(extended_y_edges) - 1), np.max(free_energy))  # 填充最大能量

        # 将原始自由能矩阵放置到新的矩阵的中心位置
        x_center_start = (len(extended_x_edges) - 1 - len(x_edges)) // 2
        y_center_start = (len(extended_y_edges) - 1 - len(y_edges)) // 2
        extended_free_energy[x_center_start:x_center_start + len(free_energy),
                                            y_center_start:y_center_start + len(free_energy[0])] = free_energy  # 将原始矩阵放在中心位置
        X, Y = np.meshgrid(extended_x_edges[:-1], extended_y_edges[:-1]) # 重新生成网格坐标用于扩展后的区域
        Z_smooth = gaussian_filter(extended_free_energy.T, sigma=1)  # 对扩展后的自由能矩阵进行高斯滤波平滑

    # 定义彩虹到白色的颜色映射
    colors = [(1, 0, 1), (0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0), (1, 0.5, 0), (1, 0, 0), (1, 1, 1)] # 定义彩虹色到白色的颜色列表
    white_to_rainbow = LinearSegmentedColormap.from_list('white_to_rainbow', colors, N=256) # 创建自定义颜色映射

    levels = np.linspace(np.min(Z_smooth), 0, 8)  # 定义等高线 levels，从最小值到 0 分成 8 级
    levels_fill = np.linspace(np.min(Z_smooth), np.max(Z_smooth), 80)  # 定义填充 levels，用于颜色填充的精细分级

    ax = subplot_ax if subplot_ax else plt.gca()  # 决定在哪个 axes 上绘图，如果 subplot_ax 为 None，则使用当前 axes

    # 绘制等高线和填充颜色
    contour = ax.contour(X, Y, Z_smooth, levels=levels, alpha=1, colors="black", linewidths=0.5, linestyles='solid') # 绘制等高线
    contourf = ax.contourf(X, Y, Z_smooth, levels=levels_fill, alpha=1, cmap=white_to_rainbow) # 绘制颜色填充
    # ax.clabel(contour, inline=True, fontsize=8, colors='black', fmt='%.1f')  # 移除等高线标签
    cbar = plt.colorbar(contourf, ax=ax, label=r"Free Energy (K$_{B}$T)")  # 添加颜色棒，显示自由能标尺，注意指定 ax
    cbar.set_ticks([])

    # 设置标题和轴标签
    if fel_type == "PCA": # PCA 类型的 FEL 图
        ax.set_xlabel(f"PC1 ({variance_explained[0]*100:.1f}%)", fontsize=10)  # X 轴标签，包含 PC1 和方差解释比例
        ax.set_ylabel(f"PC2 ({variance_explained[1]*100:.1f}%)", fontsize=10)  # Y 轴标签，包含 PC2 和方差解释比例
        ax.set_title("FEL: PC1 vs PC2", fontsize=10)  # 图标题
    elif fel_type == "RMSD-Rg": # RMSD-Rg 类型的 FEL 图
        ax.set_xlabel("RMSD (Å)", fontsize=10) # X 轴标签
        ax.set_ylabel("Radius of Gyration (Å)", fontsize=10) # Y 轴标签
        ax.set_title("FEL: RMSD vs Rg", fontsize=10) # 图标题
    ax.tick_params(axis='both', labelsize=10)  # 设置坐标轴刻度字体大小


# ---------------------------
# 修改后的 analyze_trajectory 函数
# ---------------------------


def analyze_trajectory(pro_nucle_gro_fit, pro_nucle_xtc_fit, output_dir, output_filename_prefix, image_output_dir):
    """
    分析分子动力学轨迹，计算 RMSD, RMSF 和能量 Landscape (PCA & RMSD-Rg)。

    参数:
    pro_nucle_gro_fit (str): fit 后的 GRO 文件路径。
    pro_nucle_xtc_fit (str): fit 后的 XTC 文件路径。
    output_dir (str): 输出文件保存的目录。
    output_filename_prefix (str): 输出文件名的前缀。
    image_output_dir (str): 图片输出的目录。
    """
    print(f"开始分析轨迹文件: GRO={pro_nucle_gro_fit}, XTC={pro_nucle_xtc_fit}")

    # 加载轨迹和参考结构
    try:
        u = load_trajectory(pro_nucle_gro_fit, pro_nucle_xtc_fit) # 使用新的 load_trajectory 函数
        ref = load_trajectory(pro_nucle_gro_fit, pro_nucle_gro_fit)  # 使用初始结构作为参考，并使用 load_trajectory
    except Exception as e:
        print(f"加载轨迹文件失败: {e}")
        return

    nucleic_sel = "nucleic" # 核酸原子选择语句，与notebook保持一致
    nucleic = u.select_atoms(nucleic_sel)
    ref_nucleic = ref.select_atoms(nucleic_sel)

    # 创建一个包含 1 行 4 列子图的 figure 和 axes 对象
    fig, axes = plt.subplots(1, 4, figsize=figure_size, facecolor='white')  # 使用全局变量 figure_size，设置背景色为白色

    # 1. 计算 RMSD 随时间的变化并绘图 (保持原有RMSD计算和绘图)
    print("开始计算 RMSD...")
    mobile_atoms_rmsd = u.select_atoms(nucleic_sel)
    reference_atoms_rmsd = ref.select_atoms(nucleic_sel)

    rmsd_analysis = RMSD(mobile_atoms_rmsd, reference_atoms_rmsd,
                          center=True,
                          superposition=True)
    try:
        rmsd_analysis.run()
        rmsd_values = rmsd_analysis.rmsd[:, 2]
        time_steps = range(len(rmsd_values))
        time_ns = [step * TIME_PER_STEP_NS * FRAME_SAMPLING_INTERVAL for step in time_steps] # 将 time_steps 转换为 time_ns

        axes[0].plot(time_ns, rmsd_values)  # 使用 time_ns 作为横轴数据
        axes[0].set_xlabel('Time (ns)')      # 更新 x 轴标签为 'Time (ns)'
        axes[0].set_ylabel('RMSD (Å)')
        axes[0].set_title('RMSD of Nucleic Over Time')
        print(f"RMSD 随时间变化图已准备")
    except Exception as e:
        print(f"RMSD 计算失败: {e}")


    # 2. 计算 RMSF 并绘图 (保持原有RMSF计算和绘图)
    print("开始计算 RMSF...")
    try:
        aligner_rmsf = align.AlignTraj(u, ref, select='protein', in_memory=True).run()  # 对齐用于RMSF计算
        rmsf_analyzer = RMSF(nucleic).run()

        residues = nucleic.residues.resids
        res_rmsf = []
        for res in nucleic.residues:
            atom_indices = res.atoms.ix - nucleic[0].ix
            res_rmsf.append(rmsf_analyzer.results.rmsf[atom_indices].mean())

        axes[1].plot(residues, res_rmsf, 'o-')  # 使用 axes[1] 绘制
        axes[1].set_xlabel('Residue ID')
        axes[1].set_ylabel('RMSF (Å)')
        axes[1].set_title('RMSF per Residue (Nucleic)')
        print(f"RMSF per Residue 图已准备")
    except Exception as e:
        print(f"RMSF 计算失败: {e}")


    # ---------------------------
    #  替换原有能量 Landscape 部分 (使用新的函数)
    # ---------------------------
    # 3. 基于PCA的FEL (PC1 vs PC2)  (替换原有基于RMSD的能量 Landscape)
    if run_pca_analysis:
        print("\nPerforming PCA analysis and FEL calculation (sklearn)...")
        try:
            pca_sklearn_obj, pc1_sklearn, pc2_sklearn, variance_ratio_sklearn = calculate_pca_sklearn(u, nucleic_sel, n_components=2, already_aligned=False)
            free_energy_pca, xedge_pca, yedge_pca = calculate_free_energy(pc1_sklearn, pc2_sklearn, bins_x=hist_bins_pca, bins_y=hist_bins_pca) # PCA 使用 hist_bins_pca
            plot_fel(xedge_pca, yedge_pca, free_energy_pca, "PCA", variance_explained=variance_ratio_sklearn[:2], subplot_ax=axes[2], add_padding=True) # 使用 axes[2] 绘制, 并传递子图 axes，添加留白
            print(f"PCA based FEL 图已准备")

            # 保存 PCA FEL 数据 (可选)
            output_data_filename_pca_el = os.path.join(output_dir, f"{output_filename_prefix}_pca_energy_landscape_data.txt")
            # 将X_pca, Y_pca, free_energy_pca.T  flatten 后保存为三列数据
            X_pca, Y_pca = np.meshgrid(xedge_pca[:-1], yedge_pca[:-1]) # 确保 X_pca, Y_pca 已定义
            np.savetxt(output_data_filename_pca_el, np.column_stack([X_pca.flatten(), Y_pca.flatten(), free_energy_pca.T.flatten()]),
                        header='PC1  PC2  Free Energy (K$_{B}$T)', fmt='%10.5f')
            print(f"基于 PCA 的能量 landscape 数据已保存到文件: {output_data_filename_pca_el}")


        except Exception as e:
            print(f"PCA 分析和 FEL 计算失败: {e}")
    else:
        print("跳过 PCA 分析.")


    # 4. 基于RMSD和Rg的FEL  (替换原有基于Rg的能量 Landscape)
    print("\nCalculating RMSD and Rg based FEL...")
    try:
        rmsd_values_rg, rgyr_values = calculate_rmsd_rg(u, nucleic_sel, ref_universe=ref) # 使用新的 calculate_rmsd_rg 函数
        free_energy_rg, xedge_rg, yedge_rg = calculate_free_energy(rmsd_values_rg, rgyr_values, bins_x=hist_bins_rmsd, bins_y=hist_bins_rg) # RMSD-Rg 使用 hist_bins_rmsd 和 hist_bins_rg
        plot_fel(xedge_rg, yedge_rg, free_energy_rg, "RMSD-Rg", subplot_ax=axes[3], add_padding=False) # 使用 axes[3] 绘制, 并传递子图 axes，不加留白
        print(f"RMSD vs Rg based FEL 图已准备")

        # 保存 RMSD-Rg FEL 数据 (可选)
        output_data_filename_rg_el = os.path.join(output_dir, f"{output_filename_prefix}_rmsd_rg_energy_landscape_data.txt")
        # 将X_rg, Y_rg, free_energy_rg.T flatten 后保存为三列数据
        X_rg, Y_rg = np.meshgrid(xedge_rg[:-1], yedge_rg[:-1]) # 确保 X_rg, Y_rg 已定义
        np.savetxt(output_data_filename_rg_el, np.column_stack([X_rg.flatten(), Y_rg.flatten(), free_energy_rg.T.flatten()]),
                    header='RMSD (Å)  Rg (Å)  Free Energy (K$_{B}$T)', fmt='%10.5f')
        print(f"基于 RMSD-Rg 的能量 landscape 数据已保存到文件: {output_data_filename_rg_el}")


    except Exception as e:
        print(f"RMSD 和 Rg based FEL 计算失败: {e}")

    plt.tight_layout()  # 自动调整子图布局，避免重叠
    output_plot_filename_all = os.path.join(output_dir, f"{output_filename_prefix}_all_plots.png")
    fig.savefig(output_plot_filename_all)  # 保存包含所有子图的 figure
    print(f"所有分析图形已保存为: {output_plot_filename_all}")


    # 将生成的 all_plots 图片复制到指定的图片输出文件夹 (保持原有图片复制逻辑)
    try:
        md_dirname = os.path.basename(output_dir) # 获取 md_dir 的目录名
        image_output_path_all_plots = os.path.join(image_output_dir, f"{md_dirname}_all_plots.png") # 使用 md_dirname 作为新的文件名
        fig.savefig(image_output_path_all_plots) # 保存 all plots 图
        print(f"分析图片已复制到图片输出文件夹: {image_output_dir}")
    except Exception as e:
        print(f"复制分析图片到输出文件夹失败: {e}")


    # plt.show()  #  统一显示所有图形 # 注释掉，不需要显示图形

    print(f"轨迹分析完成。结果文件已保存到: {output_dir}, 图片已保存到: {image_output_dir}")



def process_trajectory(md_dir, frame_sampling_interval):
    """
    处理分子动力学轨迹，提取蛋白质和核酸，并进行fit和rotate。抽帧间隔可变。

    参数:
    md_dir (str): 分子动力学模拟目录的路径。
    frame_sampling_interval (int): 轨迹抽帧间隔，每 frame_sampling_interval 帧抽一帧。
    """
    print(f"开始处理轨迹文件，目录: {md_dir}, 抽帧间隔: {frame_sampling_interval}")

    input_PDB = os.path.join(md_dir, "outsys.pdb")
    input_gro = os.path.join(md_dir, "outsys_new.gro")
    input_xtc = os.path.join(md_dir, "output.xtc")
    pro_nucle_gro = os.path.join(md_dir, f"PD_sample{frame_sampling_interval}.gro") # 文件名包含抽帧信息
    pro_nucle_xtc = os.path.join(md_dir, f"PD_sample{frame_sampling_interval}.xtc") # 文件名包含抽帧信息
    pro_nucle_gro_fit = os.path.join(md_dir, f"PD_fit_sample{frame_sampling_interval}.gro") # 文件名包含抽帧信息
    pro_nucle_xtc_fit = os.path.join(md_dir, f"PD_fit_sample{frame_sampling_interval}.xtc") # 文件名包含抽帧信息

    # 1. PDB to GRO 转换 (如果 gro 文件不存在)
    if not os.path.exists(input_gro):
        print(f"将PDB文件转换为GRO文件: {input_PDB} -> {input_gro}")
        try:
            u_pdb = mda.Universe(input_PDB, format="pdb")
            u_pdb.atoms.write(input_gro)
        except Exception as e:
            print(f"PDB to GRO 转换失败: {e}")
            return

    # 2. 读取 gro 和 xtc 文件
    try:
        u = mda.Universe(input_gro, input_xtc)
    except Exception as e:
        print(f"读取GRO和XTC文件失败: {e}")
        return

    # 3. 选择蛋白质和核酸并写入新的 gro 和 xtc (按指定间隔抽帧)
    protein_and_nucleic = u.select_atoms('protein or nucleic')
    protein_and_nucleic.dimensions = u.dimensions # 保留盒子尺寸信息
    try:
        protein_and_nucleic.write(pro_nucle_gro) # 保存为新的 gro
        with mda.Writer(pro_nucle_xtc, protein_and_nucleic.atoms.n_atoms) as W: # 保存为新的 xtc (按指定间隔抽帧)
            for ts in u.trajectory[::frame_sampling_interval]: # 使用可变抽帧间隔
                W.write(protein_and_nucleic.atoms)
        print(f"已保存蛋白质和核酸的 gro 和 xtc 文件: {pro_nucle_gro}, {pro_nucle_xtc}")
    except Exception as e:
        print(f"保存蛋白质和核酸的 gro/xtc 文件失败: {e}")
        return


    # 4. Fit and Rotate 轨迹处理
    try:
        u_fit = mda.Universe(pro_nucle_gro, pro_nucle_xtc) # 使用新的 gro 和 xtc
        protein = u_fit.select_atoms('protein')
        dna = u_fit.select_atoms('nucleic')
        complex_fit = u_fit.select_atoms('protein or nucleic')
        u_fit.atoms.guess_bonds()
        ref_u = u_fit.copy()
        reference = ref_u.select_atoms("protein")

        workflow = (
            unwrap(complex_fit),
            center_in_box(dna, center='mass'),
            wrap(complex_fit, compound='fragments'),
            fit_rot_trans(protein, reference)
        )
        u_fit.trajectory.add_transformations(*workflow)

        u_fit.atoms.write(pro_nucle_gro_fit) # 保存 fit 后的 gro
        with mda.Writer(pro_nucle_xtc_fit, u_fit.atoms.n_atoms) as W: # 保存 fit 后的 xtc (所有帧)
            for ts in u_fit.trajectory:
                W.write(u_fit.atoms)
        print(f"已保存 fit & rotated 的 gro 和 xtc 文件: {pro_nucle_gro_fit}, {pro_nucle_xtc_fit}")

    except Exception as e:
        print(f"轨迹 Fit & Rotate 处理失败: {e}")
        return

    print(f"轨迹文件处理完成，目录: {md_dir}, 抽帧间隔: {frame_sampling_interval}")



def main():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description="分子动力学轨迹分析脚本，可指定父目录和图片输出基础目录。")

    # 添加命令行参数
    parser.add_argument('--md_parent_dir', type=str, required=True,
                        help='分子动力学模拟目录的父目录路径')
    parser.add_argument('--image_output_base_dir', type=str, required=True,
                        help='图片输出文件夹的基础目录路径')

    # 解析命令行参数
    args = parser.parse_args()

    # 从命令行参数获取 md_parent_dir 和 image_output_base_dir
    md_parent_dir = args.md_parent_dir
    image_output_base_dir = args.image_output_base_dir
    image_output_dir = os.path.join(image_output_base_dir, IMAGE_OUTPUT_DIR_NAME) # 图片输出文件夹名称，基于传入的基础目录

    # --- 需要处理的 md_dir 列表 ---
    # ... (目录列表获取部分，保持不变) ...
    md_dirs = []
    parent_dir_contents = os.listdir(md_parent_dir)

    for item_name in parent_dir_contents:
        item_path = os.path.join(md_parent_dir, item_name)
        if os.path.isdir(item_path):
            if item_name.startswith('iter2_'):
                md_dirs.append(item_path)
    # --- md_dir 列表配置结束 ---


    if not os.path.exists(image_output_dir):
        os.makedirs(image_output_dir)  # 如果图片输出文件夹不存在，则创建

    for md_dir in md_dirs:
        if not os.path.isdir(md_dir): # 确保是目录
            print(f"跳过非目录路径: {md_dir}")
            continue

        print(f"开始处理目录: {md_dir}")
        simulation_state_pkl = os.path.join(md_dir, "simulation_state.pkl") # 检查 simulation_state.pkl 是否存在
        pro_nucle_gro_fit = os.path.join(md_dir, f"PD_fit_sample{FRAME_SAMPLING_INTERVAL}.gro") # 检查对应抽帧间隔的 fit 后 gro 文件是否存在
        pro_nucle_xtc_fit = os.path.join(md_dir, f"PD_fit_sample{FRAME_SAMPLING_INTERVAL}.xtc") # 检查对应抽帧间隔的 fit 后 xtc 文件是否存在
        output_dir = md_dir # 输出文件路径为 md_dir 自身

        if os.path.exists(pro_nucle_gro_fit) and os.path.exists(pro_nucle_xtc_fit): # 如果 fit 后的 gro 和 xtc 文件已存在，则直接执行分析 (检查对应抽帧间隔的文件)
            print(f"发现已存在的 fit 后文件 (抽帧间隔为 {FRAME_SAMPLING_INTERVAL}): {pro_nucle_gro_fit}, {pro_nucle_xtc_fit}。跳过轨迹处理，直接进行分析。")
            analyze_trajectory(pro_nucle_gro_fit, pro_nucle_xtc_fit, output_dir, OUTPUT_FILENAME_PREFIX, image_output_dir)
        elif os.path.exists(simulation_state_pkl): # 如果 simulation_state.pkl 存在，但 fit 后文件不存在，则执行轨迹处理和分析
            print(f"发现 simulation_state.pkl 文件: {simulation_state_pkl}。执行轨迹处理和分析，抽帧间隔为 {FRAME_SAMPLING_INTERVAL}。")
            process_trajectory(md_dir, FRAME_SAMPLING_INTERVAL) # 执行轨迹处理，传入抽帧间隔参数
            analyze_trajectory(pro_nucle_gro_fit, pro_nucle_xtc_fit, output_dir, OUTPUT_FILENAME_PREFIX, image_output_dir) # 执行分析
        else: # 如果 simulation_state.pkl 和 fit 后文件都不存在，则跳过该目录
            print(f"在目录 {md_dir} 中未发现 simulation_state.pkl 或已处理的文件，跳过该目录。")
        print(f"目录 {md_dir} 处理完成。\n")

    print("所有目录处理完成。")


if __name__ == "__main__":
    main()
