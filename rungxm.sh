#!/bin/bash

# GROMACS 能量最小化、平衡化 (NVT/NPT) 和分子动力学模拟脚本
#
# 修改点：
# 1. 直接创建/覆盖 mdp 文件，不进行存在性检查。
# 2. 增加了 NVT 和 NPT 平衡化步骤。
# 3. 增加了生产分子动力学模拟步骤。
# 4. 使用变量定义模拟参数和文件前缀。
# 5. 确保每个步骤的输入结构文件是上一步的输出。
# 6. 在生产模拟中禁用位置限制。
# 7. 提取最终生产模拟后的蛋白质和 DNA 结构。

# -----------------------------------------------------------------------------
# 用户可修改参数
# -----------------------------------------------------------------------------

# 定义输入PDB文件名
INPUT_PDB="3hxq_complex.pdb"

# 定义输出文件前缀
# 能量最小化: em
# NVT 平衡化: nvt
# NPT 平衡化: npt
# 生产模拟: md
OUTPUT_PREFIX_EM="em"
OUTPUT_PREFIX_NVT="nvt"
OUTPUT_PREFIX_NPT="npt"
OUTPUT_PREFIX_MD="md"

# 定义力场 (根据你的体系选择合适的力场，例如amber99sb-ildn, charmm36-jul2022等)
# 你可能需要根据你的GROMACS安装找到力场的实际名称
FORCEFIELD="amber99sb-ildn"

# 定义水模型 (例如 spc, tip3p, tip4p等)
WATER_MODEL="tip3p"

# 定义盒子类型 (例如 cubic, dodecahedron, triclinic)
BOX_TYPE="dodecahedron"

# 定义盒子边界与体系的距离 (nm)
BOX_DISTANCE=1.0

# 定义模拟温度 (Kelvin)
TEMPERATURE=300

# 定义模拟步长 (ps)
# 2 fs = 0.002 ps
DT=0.002

# 定义能量最小化最大步数
EM_STEPS=50000

# 定义 NVT 平衡化步数 (例如 100 ps @ 2 fs/step = 50000 steps)
NVT_STEPS=50000

# 定义 NPT 平衡化步数 (例如 100 ps @ 2 fs/step = 50000 steps)
NPT_STEPS=50000

# 定义生产模拟步数 (例如 1 ns @ 2 fs/step = 500000 steps)
# !!! 重要：这通常需要非常大的值 (例如 10^6 到 10^9 步) !!!
MD_STEPS=500000

# 定义能量最小化参数文件 (.mdp) 的名称
MDP_FILE_EM="minim.mdp"
MDP_FILE_IONS="ions.mdp"

# 定义 NVT 平衡化参数文件 (.mdp) 的名称
MDP_FILE_NVT="nvt.mdp"

# 定义 NPT 平衡化参数文件 (.mdp) 的名称
MDP_FILE_NPT="npt.mdp"

# 定义生产模拟参数文件 (.mdp) 的名称
MDP_FILE_MD="md.mdp"

# 定义最终输出的蛋白质+DNA PDB 文件名
# 将从最终生产模拟的结构中提取
PROTEIN_DNA_PDB="${OUTPUT_PREFIX_MD}_protein_dna.pdb"


# -----------------------------------------------------------------------------
# 脚本开始执行
# -----------------------------------------------------------------------------
set -e # 如果任何命令失败，立即退出脚本

echo "--- GROMACS MD 模拟脚本 ---"
echo "输入 PDB: ${INPUT_PDB}"
echo "力场: ${FORCEFIELD}"
echo "水模型: ${WATER_MODEL}"
echo "模拟温度: ${TEMPERATURE} K"
echo "模拟步长: ${DT} ps"
echo "EM 步数: ${EM_STEPS}"
echo "NVT 步数: ${NVT_STEPS}"
echo "NPT 步数: ${NPT_STEPS}"
echo "MD 步数: ${MD_STEPS}"
echo "------------------------------------"

# -----------------------------------------------------------------------------
# 步骤 0: 创建或覆盖 mdp 文件
# -----------------------------------------------------------------------------
echo "正在创建或覆盖 mdp 文件..."

# ions.mdp (用于加离子前的 grompp)
cat << EOF > ${MDP_FILE_IONS}
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = ${EM_STEPS}         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	= Verlet    ; Buffered neighbor searching 
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF
echo "${MDP_FILE_IONS} 创建成功."

# minim.mdp (用于能量最小化)
cat << EOF > ${MDP_FILE_EM}
; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = ${EM_STEPS}         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF
echo "${MDP_FILE_EM} 创建成功."

# nvt.mdp (用于 NVT 平衡化)
cat << EOF > ${MDP_FILE_NVT}
title                   = OPLS Lysozyme NVT equilibration 
define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
dt          = ${DT}         ; Time step (ps)
nsteps      = ${NVT_STEPS}  ; Number of steps
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = ${TEMPERATURE}     ${TEMPERATURE}           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = ${TEMPERATURE}      ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
EOF
echo "${MDP_FILE_NVT} 创建成功."


# npt.mdp (用于 NPT 平衡化)
cat << EOF > ${MDP_FILE_NPT}
title                   = OPLS Lysozyme NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
dt                      = ${DT}         ; Time step (ps)
nsteps                  = ${NPT_STEPS}  ; Number of steps
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes       ; Restarting after NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = ${TEMPERATURE}     ${TEMPERATURE}           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no        ; Velocity generation is off 
EOF
echo "${MDP_FILE_NPT} 创建成功."


# md.mdp (用于生产模拟)
cat << EOF > ${MDP_FILE_MD}
title                   = OPLS Lysozyme NPT equilibration 
; Run parameters
integrator              = md        ; leap-frog integrator
dt          = ${DT}         ; Time step (ps)
nsteps      = ${MD_STEPS}   ; Number of steps
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System    ; save the whole system
; Bond parameters
continuation            = yes       ; Restarting after NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = ${TEMPERATURE}     ${TEMPERATURE}           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel                 = no        ; Velocity generation is off 
EOF
echo "${MDP_FILE_MD} 创建成功."

echo "mdp 文件创建完毕."

# -----------------------------------------------------------------------------
# 步骤 1: 生成拓扑文件和处理过的结构文件
# gmx pdb2gmx: 从PDB文件生成拓扑文件 (.top) 和位置限制文件 (.itp)
# -f: 输入PDB文件
# -o: 输出GROMACS格式的结构文件 (.gro)
# -p: 输出拓扑文件 (.top)
# -ff: 选择力场
# -water: 选择水模型
# -ignh: 忽略PDB文件中的氢原子，GROMACS会根据力场重新添加
# -----------------------------------------------------------------------------
echo "正在生成拓扑文件和处理过的结构文件..."
gmx pdb2gmx -f ${INPUT_PDB} -o processed.gro -p topol.top -ff ${FORCEFIELD} -water ${WATER_MODEL} -ignh
echo "拓扑文件和处理过的结构文件生成成功：processed.gro, topol.top"

# -----------------------------------------------------------------------------
# 步骤 2: 定义模拟盒子
# gmx editconf: 定义模拟盒子的大小和形状
# -f: 输入GROMACS格式的结构文件
# -o: 输出带盒子的GROMACS格式的结构文件
# -c: 将体系放置在盒子的中心
# -d: 设置体系边界到盒子边缘的距离
# -bt: 定义盒子类型
# -----------------------------------------------------------------------------
echo "正在定义模拟盒子..."
gmx editconf -f processed.gro -o box.gro -c -d ${BOX_DISTANCE} -bt ${BOX_TYPE}
echo "模拟盒子定义成功：box.gro"

# -----------------------------------------------------------------------------
# 步骤 3: 加水溶剂化体系
# gmx solvate: 向盒子中添加溶剂分子
# -cp: 输入中心放置并定义好盒子的结构文件
# -cs: 选择溶剂构型文件 (通常与选择的水模型对应，例如 spc216.gro)
# -o: 输出溶剂化后的结构文件
# -p: 体系的拓扑文件 (需要更新以包含水分子)
# -----------------------------------------------------------------------------
echo "正在向体系中加水..."
# 注意：cs 后面的文件名取决于你的GROMACS安装和选择的水模型
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top
echo "体系溶剂化成功：solvated.gro"

# -----------------------------------------------------------------------------
# 步骤 4: 添加离子以中和体系电荷
# 首先需要一个运行 grompp 的步骤来生成一个用于 genion 的 .tpr 文件
# gmx grompp: 预处理，将结构文件、拓扑文件和mdp文件组合成一个二进制运行文件 (.tpr)
# -f: 输入mdp文件 (这里使用上面创建的 ions.mdp)
# -c: 输入溶剂化后的结构文件
# -p: 输入拓扑文件
# -o: 输出.tpr文件
# -----------------------------------------------------------------------------
echo "正在准备添加离子..."
gmx grompp -f ${MDP_FILE_IONS} -c solvated.gro -p topol.top -o ions.tpr -maxwarn 5
echo ".tpr 文件生成成功：ions.tpr"

# gmx genion: 向体系中添加离子
# -s: 输入.tpr文件
# -o: 输出加离子后的结构文件
# -p: 体系的拓扑文件 (需要更新以包含离子)
# -pname: 正离子名称 (根据力场和水模型确定，例如 NA, K)
# -nname: 负离子名称 (根据力场和水模型确定，例如 CL)
# -neutral: 中和体系的总电荷
# !!! 重要提示：此处输入的组编号 '14' 是溶剂组 (SOL) 的常见编号猜测。
# !!! 您必须根据 'gmx genion -s ions.tpr' 运行时实际列出的组列表来确认并修改此数字！
# !!! 运行 'gmx genion -s ions.tpr' 交互式地运行一次，查看组编号，然后将 '14' 替换为正确的编号。
# !!! 否则，此步骤可能会失败或将离子添加到错误的组中。
echo "正在添加离子以中和体系电荷..."
gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral << EOF
14
EOF
echo "离子添加成功：ionized.gro"


# -----------------------------------------------------------------------------
# 步骤 5: 能量最小化
# -----------------------------------------------------------------------------
echo "--- 开始能量最小化 ---"
# gmx grompp: 预处理能量最小化步骤
# -f: 输入能量最小化mdp文件
# -c: 输入加离子后的结构文件
# -p: 输入拓扑文件
# -o: 输出能量最小化运行文件 (.tpr)
# -----------------------------------------------------------------------------
echo "正在预处理能量最小化步骤..."
gmx grompp -f ${MDP_FILE_EM} -c ionized.gro -p topol.top -o ${OUTPUT_PREFIX_EM}.tpr
echo "能量最小化 .tpr 文件生成成功：${OUTPUT_PREFIX_EM}.tpr"

# -----------------------------------------------------------------------------
# gmx mdrun: 运行能量最小化
# -v: 详细输出
# -deffnm: 定义输出文件名的前缀
# -----------------------------------------------------------------------------
echo "正在运行能量最小化..."
gmx mdrun -v -deffnm ${OUTPUT_PREFIX_EM}
echo "能量最小化完成！输出文件前缀：${OUTPUT_PREFIX_EM}"


# -----------------------------------------------------------------------------
# 步骤 6: NVT 平衡化
# -----------------------------------------------------------------------------
echo "--- 开始 NVT 平衡化 ---"
# gmx grompp: 预处理 NVT 步骤
# -f: 输入 NVT mdp 文件
# -c: 输入能量最小化后的结构文件
# -p: 输入拓扑文件
# -r: 输入参考结构文件用于位置限制 (通常是能量最小化后的结构)
# -o: 输出 NVT 运行文件 (.tpr)
# -----------------------------------------------------------------------------
echo "正在预处理 NVT 平衡化步骤..."
gmx grompp -f ${MDP_FILE_NVT} -c ${OUTPUT_PREFIX_EM}.gro -r ${OUTPUT_PREFIX_EM}.gro -p topol.top -o ${OUTPUT_PREFIX_NVT}.tpr -maxwarn 5
echo "NVT 平衡化 .tpr 文件生成成功：${OUTPUT_PREFIX_NVT}.tpr"

# -----------------------------------------------------------------------------
# gmx mdrun: 运行 NVT 平衡化
# -v: 详细输出
# -deffnm: 定义输出文件名的前缀
# -----------------------------------------------------------------------------
echo "正在运行 NVT 平衡化..."
# 使用 -nt 或 -ntmpi, -ntomp 控制并行运行，例如 gmx mdrun -nt 8 ...
gmx mdrun -v -deffnm ${OUTPUT_PREFIX_NVT}
echo "NVT 平衡化完成！输出文件前缀：${OUTPUT_PREFIX_NVT}"


# -----------------------------------------------------------------------------
# 步骤 7: NPT 平衡化
# -----------------------------------------------------------------------------
echo "--- 开始 NPT 平衡化 ---"
# gmx grompp: 预处理 NPT 步骤
# -f: 输入 NPT mdp 文件
# -c: 输入 NVT 结束时的结构文件 (mdrun 会输出 ${OUTPUT_PREFIX_NVT}.gro)
# -p: 输入拓扑文件
# -r: 输入参考结构文件用于位置限制 (继续使用能量最小化后的结构)
# -o: 输出 NPT 运行文件 (.tpr)
# -----------------------------------------------------------------------------
echo "正在预处理 NPT 平衡化步骤..."
gmx grompp -f ${MDP_FILE_NPT} -c ${OUTPUT_PREFIX_NVT}.gro -r ${OUTPUT_PREFIX_EM}.gro -p topol.top -o ${OUTPUT_PREFIX_NPT}.tpr -maxwarn 5
echo "NPT 平衡化 .tpr 文件生成成功：${OUTPUT_PREFIX_NPT}.tpr"

# -----------------------------------------------------------------------------
# gmx mdrun: 运行 NPT 平衡化
# -v: 详细输出
# -deffnm: 定义输出文件名的前缀
# -----------------------------------------------------------------------------
echo "正在运行 NPT 平衡化..."
# 使用 -nt 或 -ntmpi, -ntomp 控制并行运行
gmx mdrun -v -deffnm ${OUTPUT_PREFIX_NPT}
echo "NPT 平衡化完成！输出文件前缀：${OUTPUT_PREFIX_NPT}"


# -----------------------------------------------------------------------------
# 步骤 8: 生产分子动力学模拟
# -----------------------------------------------------------------------------
echo "--- 开始生产模拟 ---"
# gmx grompp: 预处理生产模拟步骤
# -f: 输入生产模拟 mdp 文件 (此时 mdp 文件中应已禁用位置限制)
# -c: 输入 NPT 结束时的结构文件 (mdrun 会输出 ${OUTPUT_PREFIX_NPT}.gro)
# -p: 输入拓扑文件
# -o: 输出生产模拟运行文件 (.tpr)
# !!! 注意：此处不再使用 -r 选项，且 md.mdp 中应已禁用 -DPOSRES !!!
# -----------------------------------------------------------------------------
echo "正在预处理生产模拟步骤..."
gmx grompp -f ${MDP_FILE_MD} -c ${OUTPUT_PREFIX_NPT}.gro -p topol.top -o ${OUTPUT_PREFIX_MD}.tpr -maxwarn 5
echo "生产模拟 .tpr 文件生成成功：${OUTPUT_PREFIX_MD}.tpr"

# -----------------------------------------------------------------------------
# gmx mdrun: 运行生产模拟
# -v: 详细输出
# -deffnm: 定义输出文件名的前缀
# -----------------------------------------------------------------------------
echo "正在运行生产模拟..."
# 这是计算量最大的步骤，通常会使用并行计算。
# 示例：MPI并行 (如果你使用MPI版本GROMACS)
# mpirun -np 8 gmx mdrun -v -deffnm ${OUTPUT_PREFIX_MD} -ntomp 4
# 示例：OpenMP/线程并行
# gmx mdrun -v -deffnm ${OUTPUT_PREFIX_MD} -nt 8
# 或者让 GROMACS 自动选择最佳并行方式
gmx mdrun -v -deffnm ${OUTPUT_PREFIX_MD}
echo "生产模拟完成！输出文件前缀：${OUTPUT_PREFIX_MD}"


# -----------------------------------------------------------------------------
# 步骤 9: 从最终生产结构中提取蛋白质和 DNA 并保存为 PDB 文件
# gmx trjconv: 转换/处理轨迹和结构文件
# -f: 输入结构文件 (最终生产模拟后的结构)
# -s: 输入运行文件 (.tpr) 用于提供拓扑和组定义
# -o: 输出 PDB 文件
# -dump 0: 提取第一帧 (因为 .gro 文件只有一帧)
# 您将被提示选择要写入输出文件的组。
# -----------------------------------------------------------------------------
echo "正在从最终生产结构中提取蛋白质和 DNA..."
# !!! 重要提示：此处输入的组编号 '1 12' 是一个占位符，代表蛋白质和 DNA 的组。
# !!! 您必须根据 'gmx trjconv -f ${OUTPUT_PREFIX_MD}.gro -s ${OUTPUT_PREFIX_MD}.tpr'
# !!! 交互式地运行一次，查看 gmx trjconv 列出的组列表，找到包含蛋白质和 DNA
# !!! 的组合组（例如 'Protein_DNA', 'Protein_other' 等）的实际编号。
# !!! 将 '1 12' 替换为该实际编号(们)。
# !!! 如果没有合适的组合组，您可能需要先使用 'gmx make_ndx -f ${OUTPUT_PREFIX_MD}.tpr -o index.ndx'
# !!! 工具手动创建一个包含蛋白质和 DNA 的自定义组，然后在 trjconv 命令中用 -n index.ndx
# !!! 指定索引文件，并选择你创建的组的编号。
gmx trjconv -f ${OUTPUT_PREFIX_MD}.gro -s ${OUTPUT_PREFIX_MD}.tpr -o ${PROTEIN_DNA_PDB} -dump 0 << EOF
1 12
EOF
echo "蛋白质和 DNA 结构已成功提取到文件：${PROTEIN_DNA_PDB}"


# -----------------------------------------------------------------------------
# 主要输出文件总结：
# ${OUTPUT_PREFIX_EM}.log: 能量最小化日志文件
# ${OUTPUT_PREFIX_NVT}.log: NVT 平衡化日志文件
# ${OUTPUT_PREFIX_NPT}.log: NPT 平衡化日志文件
# ${OUTPUT_PREFIX_MD}.log: 生产模拟日志文件
#
# ${OUTPUT_PREFIX_EM}.edr: 能量最小化能量数据
# ${OUTPUT_PREFIX_NVT}.edr: NVT 能量数据
# ${OUTPUT_PREFIX_NPT}.edr: NPT 能量数据
# ${OUTPUT_PREFIX_MD}.edr: 生产模拟能量数据
#
# ${OUTPUT_PREFIX_EM}.gro: 能量最小化后的完整体系结构
# ${OUTPUT_PREFIX_NVT}.gro: NVT 结束时的完整体系结构
# ${OUTPUT_PREFIX_NPT}.gro: NPT 结束时的完整体系结构
# ${OUTPUT_PREFIX_MD}.gro: 生产模拟结束时的完整体系结构
#
# ${OUTPUT_PREFIX_NVT}.xtc: NVT 轨迹 (压缩)
# ${OUTPUT_PREFIX_NPT}.xtc: NPT 轨迹 (压缩)
# ${OUTPUT_PREFIX_MD}.xtc: 生产模拟轨迹 (压缩)
#
# ${PROTEIN_DNA_PDB}: 最终生产模拟后的蛋白质和 DNA 结构 (PDB 格式)
#
# 您可以使用 gmx energy 查看能量/温度/压力随时间的变化，例如：
# echo "Potential" | gmx energy -f ${OUTPUT_PREFIX_EM}.edr -o em_potential.xvg
# echo "Temperature" | gmx energy -f ${OUTPUT_PREFIX_NVT}.edr -o nvt_temp.xvg
# echo "Pressure" | gmx energy -f ${OUTPUT_PREFIX_NPT}.edr -o npt_press.xvg
# 然后使用 xmgrace 或其他绘图工具查看 .xvg 文件。
#
# 使用 VMD 或 PyMOL 查看结构和轨迹：
# vmd ${PROTEIN_DNA_PDB}
# vmd ${OUTPUT_PREFIX_MD}.tpr -xtc ${OUTPUT_PREFIX_MD}.xtc
# -----------------------------------------------------------------------------

echo "--- 脚本执行完毕 ---"
echo "模拟流程完成：能量最小化 -> NVT 平衡化 -> NPT 平衡化 -> 生产模拟"
echo "最终生产模拟轨迹文件: ${OUTPUT_PREFIX_MD}.xtc"
echo "最终生产模拟结束时的结构文件 (完整体系): ${OUTPUT_PREFIX_MD}.gro"
echo "请检查日志文件 (.log) 和能量文件 (.edr) 以确认模拟的稳定性。"
