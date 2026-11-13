#!/bin/bash
# ============================================================
# SATURN training script - Mouse x Human Integration
# ============================================================
# Before running, please confirm:
#   1. conda activate saturn_env
#   2. module load cuda12.1/toolkit/12.1.1
#   3. protein embeddings have been downloaded/generated
#   4. run_saturn_mouse_human.py has finished running
# ============================================================

set -e

WORK_DIR="/projects/hmz-aging/Zebrafish_Mouse_Human_inte"
SATURN_CODE="/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN"  # SATURN repo path
SATURN_DATA="$WORK_DIR/saturn_data"

cd "$SATURN_CODE"

echo "========================================="
echo "Checking GPU"
echo "========================================="
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
nvidia-smi

echo ""
echo "========================================="
echo "Checking input files"
echo "========================================="
echo "Input CSV:"
cat "$SATURN_DATA/mouse_human_run.csv"
echo ""
echo "Cell type map:"
cat "$SATURN_DATA/mouse_human_ct_map.csv"

echo ""
echo "========================================="
echo "Step 1: Download protein embeddings (if missing)"
echo "========================================="
EMB_DIR="$SATURN_DATA/protein_embeddings"
mkdir -p "$EMB_DIR"

MOUSE_EMB="/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN/protein_embeddings/protein_embeddings_export/ESM2/mouse_embedding.torch"
HUMAN_EMB="/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN/protein_embeddings/protein_embeddings_export/ESM2/human_embedding.torch"

if [ ! -f "$MOUSE_EMB" ] || [ ! -f "$HUMAN_EMB" ]; then
    echo "Protein embeddings not found, downloading..."
    cd "$EMB_DIR"
    wget -c http://snap.stanford.edu/saturn/data/protein_embeddings.tar.gz
    tar -xzf protein_embeddings.tar.gz
    echo "Download complete; please check the files and place the ESM2 .pt files in:"
    echo "  $MOUSE_EMB"
    echo "  $HUMAN_EMB"
    echo ""
    echo "Downloaded file structure:"
    find . -name "*.pt" | head -20
    cd "$SATURN_CODE"
fi

echo ""
echo "========================================="
echo "Step 2: Run SATURN training"
echo "========================================="
# --in_data:             input CSV file
# --in_label_col:        celltype label column name
# --ref_label_col:       reference label column (can be the same as in_label_col)
# --num_macrogenes:      number of macrogenes (default 2000)
# --hv_genes:            number of highly variable genes (default 8000)
# --centroids_init_path: save/load KMeans-initialized centroids
# --embedding_model:     ESM2
# --pretrain_epochs:     pretraining epochs
# --epochs:              metric learning epochs
# --device_num:          GPU index
# --work_dir:            output directory
# --score_adata:         score during training
# --ct_map_path:         celltype mapping file

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
echo "Training complete!"
echo "========================================="
echo "Output files are in: $SATURN_DATA/saturn_results/"
ls -la "$SATURN_DATA/saturn_results/"
