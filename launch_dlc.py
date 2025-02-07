#!/usr/bin/env python3

import os
import subprocess
import argparse
import time

# 动态地将 PDB 处理任务分配到 GPU 上，并且可以指定要处理的父目录。
def find_free_gpu(gpu_status):
    """查找空闲 GPU 返回第一个空闲 GPU 的索引，如果没有则返回 None"""
    for i, status in enumerate(gpu_status):
        if status == 0:
            return i
    return None

def launch_task(folder, gpu_id, gpu_status, process_pids):
    """启动 runopenmm.py 任务"""
    pdb_file = None
    for filename in os.listdir(folder):
        if filename.endswith(".pdb"):
            pdb_file = os.path.join(folder, filename)
            break

    if pdb_file:
        print(f"发现 PDB 文件: {pdb_file}, 分配到 GPU {gpu_id}, 文件夹: {folder}")
        runopenmm_path = "/path/to/runopenmm.py"
        command = ["python", runopenmm_path, pdb_file, "--gpu", str(gpu_id)]
        print(f"执行命令: {' '.join(command)}")
        process = subprocess.Popen(command) # 异步执行命令
        process_pids.append((process, gpu_id)) # 记录进程对象和 GPU ID
        gpu_status[gpu_id] = 1 # 标记 GPU
        print(f"GPU {gpu_id} 状态更新为: {gpu_status[gpu_id]}")
    else:
        print(f"警告: 在文件夹 {folder} 中没有找到 PDB 文件")

def main():
    # num_gpus = 1  # GPU 总数量
    parent_dir = "."  # 默认父目录为当前目录

    parser = argparse.ArgumentParser(description="动态分配 PDB 任务到 GPU 脚本")
    parser.add_argument("-d", "--dir", type=str, help="父目录路径", default=".")
    parser.add_argument("-g", "--gpus", type=int, help="GPU 总数量", default=1)
    args = parser.parse_args()
    parent_dir = args.dir
    num_gpus = args.gpus
    if not os.path.isdir(parent_dir):
        print(f"错误: 指定的父目录 '{parent_dir}' 不存在")
        return

    gpu_status = [0] * num_gpus  # 初始化 GPU 状态列表，0 表示空闲，1 表示 занято
    process_pids = []  # 存储后台任务的 (进程对象, GPU ID) 列表
    folder_queue = []  # 待处理的文件夹队列

    # 构建待处理文件夹队列
    for item in os.listdir(parent_dir):
        folder_path = os.path.join(parent_dir, item)
        if os.path.isdir(folder_path) and folder_path != parent_dir:
            folder_queue.append(folder_path)

    if not folder_queue:
        print(f"在目录 '{parent_dir}' 中没有找到任何子文件夹")
        return

    print(f"待处理文件夹队列: {folder_queue}")
    print("开始动态分配任务...")

    while True:
        # 1. 检查是否有已完成的任务，释放 GPU
        i = 0
        while i < len(process_pids):
            process, gpu_id = process_pids[i]
            if process.poll() is not None: # 检查进程是否结束
                completed_pid = process.pid
                completed_gpu_id = gpu_id

                print(f"任务 PID {completed_pid} (GPU {completed_gpu_id}) 已完成")
                gpu_status[completed_gpu_id] = 0 # 标记 GPU 为空闲
                print(f"GPU {completed_gpu_id} 状态更新为: {gpu_status[completed_gpu_id]}")

                process_pids.pop(i) # 移除已完成的任务 (注意 pop(i) 会改变列表长度，所以这里保持 i 不变，下次循环检查当前索引的新元素)
                continue # 移除元素后，保持 i 不变，继续检查当前索引，避免跳过元素
            i += 1

        # 2. 分配新任务到空闲 GPU
        while True:
            free_gpu_id = find_free_gpu(gpu_status)
            if free_gpu_id is not None: # 找到空闲 GPU
                if folder_queue:
                    current_folder = folder_queue.pop(0) # 从队列头部取出一个文件夹
                    launch_task(current_folder, free_gpu_id, gpu_status, process_pids)
                else:
                    print("所有文件夹已加入任务队列，等待剩余任务完成...")
                    break # 跳出内层循环，进入等待剩余任务完成阶段
            else:
                # 没有空闲 GPU，跳出内层循环，稍后重试
                break

        # 3. 检查是否所有任务都已分配且完成
        if not folder_queue and not process_pids:
            print("所有任务执行完毕.")
            break # 跳出主循环，脚本结束

        time.sleep(10) # 稍作等待，降低 CPU 占用，并给任务完成留出时间

    print("脚本执行完成.")

if __name__ == "__main__":
    main()
