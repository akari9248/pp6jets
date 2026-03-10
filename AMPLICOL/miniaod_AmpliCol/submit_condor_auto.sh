#!/bin/bash
# 自动提交 condor 任务的脚本（无需确认，适合批量提交）
# 用法: ./submit_condor_auto.sh <Part编号> [队列数量]
# 示例: ./submit_condor_auto.sh 1       # 提交 Part1, 队列1000（默认）
#       ./submit_condor_auto.sh 2 500   # 提交 Part2, 队列500
# 批量: for i in {1..5}; do ./submit_condor_auto.sh $i; done

set -e

# 检查参数
if [ -z "$1" ]; then
    echo "用法: $0 <Part编号> [队列数量]"
    echo "示例: $0 1       # 提交 Part1, 队列1000"
    echo "      $0 2 500   # 提交 Part2, 队列500"
    echo "批量: for i in {1..5}; do $0 \$i; done"
    exit 1
fi

PART_NUM=$1
QUEUE_NUM=${2:-1000}

# 脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDOR_JDL="${SCRIPT_DIR}/condor.jdl"

# 提交任务 - 通过命令行参数传递变量，不再修改文件！
cd "${SCRIPT_DIR}"
mkdir -p log
echo "[Part${PART_NUM}] 提交任务 (队列: ${QUEUE_NUM})..."
condor_submit "${CONDOR_JDL}" PART_NUM="${PART_NUM}" QUEUE_NUM="${QUEUE_NUM}"

echo "[Part${PART_NUM}] 完成！"

