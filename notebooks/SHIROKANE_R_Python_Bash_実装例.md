# SHIROKANE R/Python/Bash 実装例

#SHIROKANE #R #Python #Bash #実装例 #スパコン #ジョブスクリプト #AGE #パイプライン #ワークフロー

このドキュメントは、SHIROKANEスーパーコンピュータでのR、Python、Bashを連携させた生物情報学解析の実装例を示します。

## 基本ジョブスクリプト例

### シンプルなR言語ジョブ

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=4G
#$ -pe def_slot 2
#$ -o logs/r_job_$JOB_ID.log
#$ -e logs/r_job_$JOB_ID.err

# 作業ディレクトリ作成
mkdir -p logs results

# モジュールロード
module load R/4.1.0

# 実行ログ
echo "Job started: $(date)"
echo "Working directory: $(pwd)"
echo "Host: $(hostname)"

# Rスクリプト実行
Rscript scripts/analysis.R \
  --input data/expression.csv \
  --output results/differential_expression.csv

echo "Job finished: $(date)"
```

### シンプルなPythonジョブ

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=8G
#$ -pe def_slot 4
#$ -o logs/py_job_$JOB_ID.log
#$ -e logs/py_job_$JOB_ID.err

# 作業ディレクトリ作成
mkdir -p logs results

# モジュールロード
module load python/3.9

# 実行ログ
echo "Job started: $(date)"
echo "Python version: $(python --version)"

# Python実行
python scripts/process_data.py \
  --input data/sequences.fastq \
  --output results/processed_sequences.csv \
  --threads $NSLOTS

echo "Job finished: $(date)"
```

## 実際のR言語スクリプト例

### 差次的発現解析スクリプト (analysis.R)

```r
#!/usr/bin/env Rscript

# ライブラリのロード
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(data.table)
  library(futile.logger)
  library(parallel)
})

# コマンドライン引数の設定
option_list <- list(
  make_option("--input", type="character", help="Input count matrix file"),
  make_option("--metadata", type="character", help="Sample metadata file"),
  make_option("--output", type="character", help="Output results file"),
  make_option("--threads", type="integer", default=1, help="Number of CPU threads")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ロギング設定
log_file <- file.path("logs", paste0("DESeq2_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
flog.appender(appender.file(log_file))
flog.threshold(INFO)

flog.info("Starting differential expression analysis")
flog.info("Parameters: %s", paste(names(opt), unlist(opt), sep="=", collapse=", "))

# スレッド数の設定
num_threads <- opt$threads
flog.info("Using %d threads", num_threads)
BiocParallel::register(BiocParallel::MulticoreParam(num_threads))

# データの読み込み
flog.info("Reading count data from: %s", opt$input)
counts <- fread(opt$input)
gene_ids <- counts$gene_id
counts$gene_id <- NULL

flog.info("Reading metadata from: %s", opt$metadata)
metadata <- fread(opt$metadata)

# DESeq2データセットの作成
flog.info("Creating DESeq2 dataset")
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts),
  colData = metadata,
  design = ~ condition
)

# 低発現遺伝子のフィルタリング
flog.info("Filtering low-expressed genes")
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
flog.info("Retained %d genes after filtering", nrow(dds))

# 差次的発現解析の実行
flog.info("Running DESeq2 analysis")
dds <- DESeq(dds)
flog.info("Extracting results")
res <- results(dds, contrast=c("condition", "treatment", "control"))

# 結果の整形と保存
flog.info("Processing and saving results")
res_df <- as.data.frame(res)
res_df$gene_id <- gene_ids[keep]
res_df <- res_df[order(res_df$padj),]

# 結果の保存
flog.info("Writing results to: %s", opt$output)
write.csv(res_df, file=opt$output, row.names=FALSE)

flog.info("Analysis completed successfully")
```

## 実際のPythonスクリプト例

### シーケンスデータ処理スクリプト (process_data.py)

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import multiprocessing as mp

# ロギング設定
def setup_logger():
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"process_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

# コマンドライン引数の設定
def parse_args():
    parser = argparse.ArgumentParser(description="Process FASTQ sequence data")
    parser.add_argument("--input", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--min-length", type=int, default=100, help="Minimum sequence length")
    parser.add_argument("--min-quality", type=float, default=20.0, help="Minimum average quality score")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads")
    return parser.parse_args()

# シーケンスの品質チェック
def check_sequence_quality(record):
    """シーケンスの品質をチェックし、フィルタリング条件を満たすかどうかを返す"""
    avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record.seq)
    return len(record.seq) >= args.min_length and avg_quality >= args.min_quality

# シーケンスの処理（並列処理用）
def process_sequence_batch(batch):
    """シーケンスバッチを処理し、結果を返す"""
    results = []
    for record in batch:
        if check_sequence_quality(record):
            gc_content = (record.seq.count('G') + record.seq.count('C')) / len(record.seq) * 100
            results.append({
                'id': record.id,
                'length': len(record.seq),
                'gc_content': gc_content,
                'avg_quality': sum(record.letter_annotations["phred_quality"]) / len(record.seq)
            })
    return results

# メイン関数
def main():
    logger = setup_logger()
    global args
    args = parse_args()
    
    logger.info(f"Starting sequence processing")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Output file: {args.output}")
    logger.info(f"Parameters: min_length={args.min_length}, min_quality={args.min_quality}, threads={args.threads}")
    
    # 入力ファイルの存在チェック
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # 出力ディレクトリの作成
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    
    # シーケンスの読み込みと処理
    logger.info("Reading sequences from FASTQ file")
    records = list(SeqIO.parse(args.input, "fastq"))
    logger.info(f"Total sequences: {len(records)}")
    
    # バッチに分割して並列処理
    batch_size = max(1, len(records) // args.threads)
    batches = [records[i:i + batch_size] for i in range(0, len(records), batch_size)]
    logger.info(f"Processing in {len(batches)} batches with {args.threads} threads")
    
    # 並列処理の実行
    with mp.Pool(processes=args.threads) as pool:
        batch_results = pool.map(process_sequence_batch, batches)
    
    # 結果の結合
    all_results = []
    for batch_result in batch_results:
        all_results.extend(batch_result)
    
    logger.info(f"Retained {len(all_results)} sequences after filtering")
    
    # 結果のDataFrameへの変換と保存
    df = pd.DataFrame(all_results)
    logger.info(f"Saving results to {args.output}")
    df.to_csv(args.output, index=False)
    
    logger.info("Processing completed successfully")

if __name__ == "__main__":
    main()
```

## 複合ワークフローの例

### 連携ワークフロージョブスクリプト

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=16G
#$ -pe def_slot 8
#$ -o logs/workflow_$JOB_ID.log
#$ -e logs/workflow_$JOB_ID.err

# 必要なディレクトリの作成
mkdir -p logs
mkdir -p data/processed
mkdir -p results/temp
mkdir -p results/final

# 作業開始ログ
echo "=========================================="
echo "Workflow started at $(date)"
echo "Running on: $(hostname)"
echo "Working directory: $(pwd)"
echo "Allocated cores: $NSLOTS"
echo "=========================================="

# モジュールのロード
module load python/3.9
module load R/4.1.0
module load samtools/1.13
module load bedtools/2.30.0

# ステップ1: Bashによるデータ前処理
echo "Step 1: Data preprocessing [$(date)]"
bash scripts/preprocess.sh \
  --input data/raw/samples.fastq \
  --output data/processed/filtered.fastq \
  --threads $NSLOTS

if [ $? -ne 0 ]; then
  echo "ERROR: Preprocessing failed!"
  exit 1
fi

# ステップ2: Pythonによる特徴抽出
echo "Step 2: Feature extraction [$(date)]"
python scripts/extract_features.py \
  --input data/processed/filtered.fastq \
  --output data/processed/features.csv \
  --model models/feature_extractor.pkl \
  --threads $((NSLOTS/2))

if [ $? -ne 0 ]; then
  echo "ERROR: Feature extraction failed!"
  exit 1
fi

# ステップ3: Rによる統計解析
echo "Step 3: Statistical analysis [$(date)]"
Rscript scripts/analyze.R \
  --input data/processed/features.csv \
  --metadata data/metadata.csv \
  --output results/temp/stats_results.rds \
  --threads $((NSLOTS/2))

if [ $? -ne 0 ]; then
  echo "ERROR: Statistical analysis failed!"
  exit 1
fi

# ステップ4: Pythonによる可視化
echo "Step 4: Visualization [$(date)]"
python scripts/visualize.py \
  --input results/temp/stats_results.rds \
  --output results/final/figures/ \
  --format pdf

if [ $? -ne 0 ]; then
  echo "ERROR: Visualization failed!"
  exit 1
fi

# ステップ5: 結果レポート生成
echo "Step 5: Report generation [$(date)]"
Rscript scripts/generate_report.R \
  --data results/temp/stats_results.rds \
  --figures results/final/figures/ \
  --output results/final/report.html

if [ $? -ne 0 ]; then
  echo "ERROR: Report generation failed!"
  exit 1
fi

# 作業完了ログ
echo "=========================================="
echo "Workflow completed successfully at $(date)"
echo "Results available in: results/final/"
echo "=========================================="
```

## マルチジョブワークフロー制御スクリプト

```bash
#!/bin/bash
# submit_workflow.sh - ワークフロー全体をSHIROKANEに投入するスクリプト

# 作業ディレクトリの設定
PROJECT_DIR=$(pwd)
SCRIPTS_DIR="${PROJECT_DIR}/job_scripts"
LOG_DIR="${PROJECT_DIR}/logs"
mkdir -p $LOG_DIR

# ワークフロー開始ログ
echo "Starting workflow submission at $(date)" | tee ${LOG_DIR}/workflow_submission.log

# ステップ1: データのダウンロード
echo "Submitting data download job..." | tee -a ${LOG_DIR}/workflow_submission.log
JOB1=$(qsub ${SCRIPTS_DIR}/01_download_data.sh)
echo "Submitted download job: $JOB1" | tee -a ${LOG_DIR}/workflow_submission.log

# ステップ2: 品質管理とフィルタリング (データダウンロード後に実行)
echo "Submitting quality control job..." | tee -a ${LOG_DIR}/workflow_submission.log
JOB2=$(qsub -hold_jid $JOB1 ${SCRIPTS_DIR}/02_quality_control.sh)
echo "Submitted QC job: $JOB2" | tee -a ${LOG_DIR}/workflow_submission.log

# ステップ3: 特徴抽出 (品質管理後に実行)
echo "Submitting feature extraction job..." | tee -a ${LOG_DIR}/workflow_submission.log
JOB3=$(qsub -hold_jid $JOB2 ${SCRIPTS_DIR}/03_extract_features.sh)
echo "Submitted feature extraction job: $JOB3" | tee -a ${LOG_DIR}/workflow_submission.log

# ステップ4: 主成分分析と統計解析 (特徴抽出後に実行)
echo "Submitting analysis job..." | tee -a ${LOG_DIR}/workflow_submission.log
JOB4=$(qsub -hold_jid $JOB3 ${SCRIPTS_DIR}/04_analyze_data.sh)
echo "Submitted analysis job: $JOB4" | tee -a ${LOG_DIR}/workflow_submission.log

# ステップ5: 可視化とレポート生成 (解析後に実行)
echo "Submitting visualization and reporting job..." | tee -a ${LOG_DIR}/workflow_submission.log
JOB5=$(qsub -hold_jid $JOB4 ${SCRIPTS_DIR}/05_visualize_report.sh)
echo "Submitted visualization job: $JOB5" | tee -a ${LOG_DIR}/workflow_submission.log

# 完了通知ジョブ (すべてのジョブ完了後に実行)
echo "Submitting completion notification job..." | tee -a ${LOG_DIR}/workflow_submission.log
JOB6=$(qsub -hold_jid $JOB5 ${SCRIPTS_DIR}/06_notify_completion.sh)
echo "Submitted notification job: $JOB6" | tee -a ${LOG_DIR}/workflow_submission.log

echo "All jobs submitted. Use 'qstat' to monitor progress." | tee -a ${LOG_DIR}/workflow_submission.log
echo "Final results will be available in 'results/final/' upon completion." | tee -a ${LOG_DIR}/workflow_submission.log
```

## 言語間連携の実装例

### R-Pythonの連携 (reticulate使用)

```r
# R-Pythonインターフェース (r_python_interface.R)
library(reticulate)
library(futile.logger)
library(optparse)

# コマンドライン引数の設定
option_list <- list(
  make_option("--input", type="character", help="Input data file"),
  make_option("--output", type="character", help="Output results file"),
  make_option("--model", type="character", help="Python model path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ロギング設定
flog.info("Starting R-Python interface script")
flog.info("Parameters: %s", paste(names(opt), unlist(opt), sep="=", collapse=", "))

# Pythonのモジュールパスを設定
flog.info("Configuring Python environment")
use_condaenv("myenv")  # conda環境を使用する場合
# または
# use_virtualenv("path/to/venv")  # venv環境を使用する場合

# Pythonモジュールのインポート
flog.info("Importing Python modules")
np <- import("numpy")
pd <- import("pandas")
pickle <- import("pickle")

# カスタムPythonモジュールのインポート
flog.info("Importing custom Python module")
model_utils <- import_from_path("model_utils", path="src/python")

# データの読み込み
flog.info("Reading data from: %s", opt$input)
data <- read.csv(opt$input)
flog.info("Data dimensions: %d rows, %d columns", nrow(data), ncol(data))

# データをPythonに渡す
flog.info("Converting data to Python format")
py_data <- r_to_py(data)

# Pythonの学習済みモデルを読み込む
flog.info("Loading Python model from: %s", opt$model)
with(pickle$load(py_open(opt$model, "rb")) %as% model, {
  # Pythonモデルを使用して予測を実行
  flog.info("Running predictions with Python model")
  predictions <- model_utils$predict(model, py_data)
})

# 結果をRのデータフレームに変換
flog.info("Converting predictions back to R format")
predictions_r <- py_to_r(predictions)

# 結果の保存
flog.info("Saving results to: %s", opt$output)
write.csv(predictions_r, file=opt$output, row.names=FALSE)

flog.info("R-Python interface script completed successfully")
```

### Python-Rの連携 (rpy2使用)

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python-R連携スクリプト (python_r_interface.py)
import os
import sys
import argparse
import logging
from datetime import datetime
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# ロギング設定
def setup_logger():
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"py_r_interface_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

# コマンドライン引数の処理
def parse_args():
    parser = argparse.ArgumentParser(description="Python-R interface for statistical analysis")
    parser.add_argument("--input", required=True, help="Input data CSV file")
    parser.add_argument("--output", required=True, help="Output results file")
    parser.add_argument("--r-script", required=True, help="R script for analysis")
    return parser.parse_args()

def main():
    # ロガーとコマンドライン引数の設定
    logger = setup_logger()
    args = parse_args()
    
    logger.info("Starting Python-R interface")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Output file: {args.output}")
    logger.info(f"R script: {args.r_script}")
    
    # R連携の初期化
    logger.info("Initializing R interface")
    pandas2ri.activate()
    
    # 必要なRパッケージのインポート
    logger.info("Importing R packages")
    base = importr('base')
    stats = importr('stats')
    
    # データの読み込み
    logger.info(f"Reading data from {args.input}")
    try:
        data = pd.read_csv(args.input)
        logger.info(f"Data dimensions: {data.shape[0]} rows, {data.shape[1]} columns")
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        sys.exit(1)
    
    # データをRに渡す
    logger.info("Converting data to R format")
    r_data = pandas2ri.py2rpy(data)
    
    # Rの環境に変数を設定
    ro.globalenv['input_data'] = r_data
    ro.globalenv['output_file'] = args.output
    
    # Rスクリプトの実行
    logger.info(f"Running R script: {args.r_script}")
    try:
        ro.r(f'source("{args.r_script}")')
        logger.info("R script executed successfully")
    except Exception as e:
        logger.error(f"Error executing R script: {e}")
        sys.exit(1)
    
    # 結果の確認
    if os.path.exists(args.output):
        logger.info(f"Results saved to {args.output}")
    else:
        logger.warning(f"Output file not found: {args.output}")
    
    logger.info("Python-R interface completed successfully")

if __name__ == "__main__":
    main()
```

## Bash実装例（前処理スクリプト）

```bash
#!/bin/bash
# preprocess.sh - FASTQ前処理スクリプト

# エラー時に停止
set -e

# デフォルト値と引数の処理
INPUT=""
OUTPUT=""
THREADS=1
MIN_LENGTH=50
MIN_QUALITY=20

# ヘルプ関数
show_help() {
    echo "Usage: $0 --input <input.fastq> --output <output.fastq> [options]"
    echo ""
    echo "Options:"
    echo "  --input       Input FASTQ file (required)"
    echo "  --output      Output FASTQ file (required)"
    echo "  --threads     Number of CPU threads (default: 1)"
    echo "  --min-length  Minimum sequence length (default: 50)"
    echo "  --min-quality Minimum quality score (default: 20)"
    echo "  --help        Show this help message"
    exit 1
}

# 引数の解析
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            INPUT="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --min-length)
            MIN_LENGTH="$2"
            shift 2
            ;;
        --min-quality)
            MIN_QUALITY="$2"
            shift 2
            ;;
        --help)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# 必須引数のチェック
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "ERROR: Input and output files are required"
    show_help
fi

# 入力ファイルの存在チェック
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

# 出力ディレクトリの作成
OUTPUT_DIR=$(dirname "$OUTPUT")
mkdir -p "$OUTPUT_DIR"

# ログディレクトリの作成
LOG_DIR="logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/preprocess_$(date +%Y%m%d_%H%M%S).log"

# ログ関数
log() {
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a "$LOG_FILE"
}

# 作業開始ログ
log "Starting FASTQ preprocessing"
log "Input file: $INPUT"
log "Output file: $OUTPUT"
log "Parameters: threads=$THREADS, min_length=$MIN_LENGTH, min_quality=$MIN_QUALITY"

# 前処理の実行
log "Step 1: Quality checking with FastQC"
module load fastqc/0.11.9
fastqc -o "$LOG_DIR" -t "$THREADS" "$INPUT"

log "Step 2: Adapter trimming with Cutadapt"
module load cutadapt/3.4
cutadapt \
    -j "$THREADS" \
    -q "$MIN_QUALITY" \
    -m "$MIN_LENGTH" \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o "${OUTPUT_DIR}/trimmed.fastq" \
    "$INPUT" \
    > "${LOG_DIR}/cutadapt_report.txt" 2>&1

log "Step 3: Quality filtering with FASTX-Toolkit"
module load fastx/0.0.14
cat "${OUTPUT_DIR}/trimmed.fastq" | \
    fastq_quality_filter \
    -q "$MIN_QUALITY" \
    -p 90 \
    -o "$OUTPUT"

# 中間ファイルの削除
log "Cleaning up temporary files"
rm "${OUTPUT_DIR}/trimmed.fastq"

# 結果の統計
INPUT_COUNT=$(grep -c "^@" "$INPUT" || true)
OUTPUT_COUNT=$(grep -c "^@" "$OUTPUT" || true)
RETAINED_PCT=$(echo "scale=2; $OUTPUT_COUNT * 100 / $INPUT_COUNT" | bc)

log "Processing completed successfully"
log "Input sequences: $INPUT_COUNT"
log "Output sequences: $OUTPUT_COUNT"
log "Retained: ${RETAINED_PCT}%"
log "Output saved to: $OUTPUT"
```

## プロンプトテンプレート

SHIROKANE上で使用するジョブスクリプトやワークフローの設計をAIアシスタントに依頼する際のテンプレートを以下に示します：

```
## SHIROKANE解析ワークフロー設計依頼

SHIROKANEスーパーコンピュータで以下の解析を行うためのワークフローを設計してください。

### 解析の概要
[解析の目的や大まかな流れを説明]

### 入力データ
- ファイル形式: [FASTQ/BAM/CSVなど]
- サイズ: [データの概算サイズ]
- サンプル数: [処理するサンプル数]

### 解析ステップ
1. [ステップ1の説明]
2. [ステップ2の説明]
3. [ステップ3の説明]
...

### 使用したい言語/ツール
- R: [バージョン、使用したいパッケージ]
- Python: [バージョン、使用したいライブラリ]
- その他ツール: [samtools, bcftoolsなど]

### リソース要件
- メモリ要件: [推定必要メモリ]
- CPUコア数: [推定必要コア数]
- 実行時間の見積: [推定実行時間]

### 出力形式
[期待される出力ファイルの形式や内容]

### 追加要件
- [エラーハンドリングの要件]
- [並列処理の要件]
- [結果の可視化要件]
- [その他の特記事項]
``` 