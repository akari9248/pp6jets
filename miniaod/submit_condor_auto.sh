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
FULLSIM_SH="${SCRIPT_DIR}/fullsim.sh"
CONDOR_JDL="${SCRIPT_DIR}/condor.jdl"

# 配置
INPUT_BASE="/eos/cms/store/group/phys_smp/ec/zhye/pp6j_15GeV"
OUTPUT_BASE="root://eoscms.cern.ch//eos/cms/store/group/phys_smp/ec/zhye/pp6j_15GeV_MINIAOD"

echo "[Part${PART_NUM}] 开始配置..."

# 修改 fullsim.sh
sed -i "s|inputfile=\"/eos/cms/store/group/phys_smp/ec/zhye/pp6j_15GeV/Part[0-9]*/chunk|inputfile=\"${INPUT_BASE}/Part${PART_NUM}/chunk|g" "${FULLSIM_SH}"

# 修改 condor.jdl
sed -i "s|output_destination = root://eoscms.cern.ch//eos/cms/store/group/phys_smp/ec/zhye/pp6j_15GeV_MINIAOD/Part[0-9]*|output_destination = ${OUTPUT_BASE}/Part${PART_NUM}|g" "${CONDOR_JDL}"
sed -i "s|^queue [0-9]*|queue ${QUEUE_NUM}|g" "${CONDOR_JDL}"

# 确保 log 目录存在
mkdir -p "${SCRIPT_DIR}/log"

# 提交任务
cd "${SCRIPT_DIR}"
echo "[Part${PART_NUM}] 提交任务 (队列: ${QUEUE_NUM})..."
condor_submit condor.jdl

echo "[Part${PART_NUM}] 完成！"

