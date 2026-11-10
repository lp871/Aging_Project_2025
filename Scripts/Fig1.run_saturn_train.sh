#!/bin/bash
# ============================================================
# SATURN 训练脚本 - Mouse x Human Integration
# ============================================================
# 使用前请确认:
#   1. conda activate saturn_env
#   2. module load cuda12.1/toolkit/12.1.1
#   3. protein embeddings 已经下载/生成
#   4. run_saturn_mouse_human.py 已运行完成
# ============================================================

set -e

WORK_DIR="/projects/hmz-aging/Zebrafish_Mouse_Human_inte"
SATURN_CODE="/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN"  # SATURN repo 路径
SATURN_DATA="$WORK_DIR/saturn_data"

cd "$SATURN_CODE"

echo "========================================="
echo "检查 GPU"
echo "========================================="
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
nvidia-smi

echo ""
echo "========================================="
echo "检查输入文件"
echo "========================================="
echo "Input CSV:"
cat "$SATURN_DATA/mouse_human_run.csv"
echo ""
echo "Cell type map:"
cat "$SATURN_DATA/mouse_human_ct_map.csv"

echo ""
echo "========================================="
echo "Step 1: 下载 protein embeddings (如不存在)"
echo "========================================="
EMB_DIR="$SATURN_DATA/protein_embeddings"
mkdir -p "$EMB_DIR"

MOUSE_EMB="/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN/protein_embeddings/protein_embeddings_export/ESM2/mouse_embedding.torch"
HUMAN_EMB="/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN/protein_embeddings/protein_embeddings_export/ESM2/human_embedding.torch"

if [ ! -f "$MOUSE_EMB" ] || [ ! -f "$HUMAN_EMB" ]; then
    echo "Protein embeddings 不存在，正在下载..."
    cd "$EMB_DIR"
    wget -c http://snap.stanford.edu/saturn/data/protein_embeddings.tar.gz
    tar -xzf protein_embeddings.tar.gz
    echo "下载完成，请检查文件并将 ESM2 的 .pt 文件放到:"
    echo "  $MOUSE_EMB"
    echo "  $HUMAN_EMB"
    echo ""
    echo "下载的文件结构:"
    find . -name "*.pt" | head -20
    cd "$SATURN_CODE"
fi

echo ""
echo "========================================="
echo "Step 2: 运行 SATURN 训练"
echo "========================================="
# --in_data:             输入 CSV 文件
# --in_label_col:        celltype 标签列名
# --ref_label_col:       参考标签列 (同 in_label_col 即可)
# --num_macrogenes:      macrogene 数量 (default 2000)
# --hv_genes:            高变基因数量 (default 8000)
# --centroids_init_path: 保存/加载 KMeans 初始化的 centroids
# --embedding_model:     ESM2
# --pretrain_epochs:     预训练 epochs
# --epochs:              metric learning epochs
# --device_num:          GPU 编号
# --work_dir:            输出目录
# --score_adata:         训练时评分
# --ct_map_path:         celltype mapping 文件

python3 train-saturn.py \
    --in_data="$SATURN_DATA/mouse_human_run.csv" \
    --in_label_col=celltype \
    --ref_label_col=celltype \
    --num_macrogenes=2000 \
    --hv_genes=8000 \
    --centroids_init_path="$SATURN_DATA/mouse_human_centroids.pkl" \
    --embedding_model=ESM2 \
    --pretrain_epochs=200 \
    --epochs=50 \
    --device_num=0 \
    --work_dir="$SATURN_DATA" \
    --score_adata \
    --ct_map_path="$SATURN_DATA/mouse_human_ct_map.csv"

echo ""
echo "========================================="
echo "训练完成！"
echo "========================================="
echo "输出文件在: $SATURN_DATA/saturn_results/"
ls -la "$SATURN_DATA/saturn_results/"
