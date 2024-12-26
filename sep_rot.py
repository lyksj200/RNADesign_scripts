from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.vectors import Vector
import numpy as np
import math
# a script for a PDB file with A and B chains, to separate and rotate B chain
 
# 定义旋转函数
def rotate_chain(chain, angle_degrees, axis='x'):
    """
    旋转链中的所有原子
    :param chain: Bio.PDB.Chain.Chain对象
    :param angle_degrees: 旋转角度（度）
    :param axis: 旋转轴，'x', 'y', 或 'z'
    """
    angle_rad = math.radians(angle_degrees)
    cos_theta = math.cos(angle_rad)
    sin_theta = math.sin(angle_rad)

    if axis == 'x':
        rotation_matrix = np.array([
            [1, 0, 0],
            [0, cos_theta, -sin_theta],
            [0, sin_theta, cos_theta]
        ])
    elif axis == 'y':
        rotation_matrix = np.array([
            [cos_theta, 0, sin_theta],
            [0, 1, 0],
            [-sin_theta, 0, cos_theta]
        ])
    elif axis == 'z':
        rotation_matrix = np.array([
            [cos_theta, -sin_theta, 0],
            [sin_theta, cos_theta, 0],
            [0, 0, 1]
        ])
    else:
        raise ValueError("轴必须是'x', 'y', 或 'z'")

    for residue in chain:
        for atom in residue:
            atom.transform(rotation_matrix, [0, 0, 0])

# 定义平移函数
def translate_chain(chain, translation_vector):
    """
    平移链中的所有原子
    :param chain: Bio.PDB.Chain.Chain对象
    :param translation_vector: 平移向量 (x, y, z)
    """
    for residue in chain:
        for atom in residue:
            atom.transform(np.identity(3), translation_vector)

# 计算质心
def calculate_centroid(chain):
    """
    计算链的质心
    :param chain: Bio.PDB.Chain.Chain对象
    :return: 质心 (x, y, z)
    """
    centroid = np.zeros(3)
    n_atoms = 0
    for residue in chain:
        for atom in residue:
            centroid += atom.coord
            n_atoms += 1
    return centroid / n_atoms

# 主函数
def separate_and_rotate_chains(pdb_file, output_file, distance=10.0, angle=45, axis='x'):
    """
    将PDB文件中的两条链分开并旋转其中一条链
    :param pdb_file: 输入PDB文件路径
    :param output_file: 输出PDB文件路径
    :param distance: 两条链之间的距离 (Å)
    :param angle: 旋转角度 (度)
    :param axis: 旋转轴，'x', 'y', 或 'z'
    """
    parser = PDBParser()
    structure = parser.get_structure('model', pdb_file)

    # 假设有两条链，A和B
    chain_A = structure[0]['A']
    chain_B = structure[0]['B']

    # 计算两条链的质心
    centroid_A = calculate_centroid(chain_A)
    centroid_B = calculate_centroid(chain_B)

    # 计算平移向量
    translation_vector = centroid_B - centroid_A
    translation_vector = translation_vector / np.linalg.norm(translation_vector) * distance

    # 平移链B
    translate_chain(chain_B, translation_vector)

    # 旋转链B
    rotate_chain(chain_B, angle, axis)

    # 保存修改后的结构
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

# 使用示例
separate_and_rotate_chains('input.pdb', 'output.pdb', distance=15.0, angle=45, axis='x')