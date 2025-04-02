# SHIROKANE スーパーコンピュータでのR/Python/Bashによる生物情報学解析ガイド

#SHIROKANE #スパコン #開発手法 #テスト #ガイドライン #rules_for_ai #composer #cursor #実践 #R #Python #Bash #renv #conda #パイプライン #ワークフロー #FAIR #git #バージョン管理 #ジョブスクリプト #AGE #TDD #RDD #要件定義

> [!NOTE]
> このガイドはSHIROKANE（東京大学ヒューマンゲノムセンター）スーパーコンピュータ上でR、Python、Bashを用いた生物情報学解析を行うためのルールをまとめたものです。
> 関連ガイド: [[../../SHIROKANE_tip/SHIROKANE_インデックス|SHIROKANEインデックス]]

## 1. 基本原則

- 解析の再現性を最優先し、すべての手順を明示的に記録する
- ジョブスクリプトを使用して計算リソースを効率的に活用する
- コードとデータを適切に分離し、バージョン管理する
- テスト駆動開発（TDD）の原則を取り入れ、解析の堅牢性を確保する
- FAIR原則（Findable, Accessible, Interoperable, Reusable）に準拠する
- 異なる言語間の連携を明確に定義し、効率的なワークフローを構築する
- ユーザーとの対話を通じて要件を明確に定義し、仕様書を作成してから開発を進める（要件駆動開発）

## 1.1 要件定義と開発プロセス

解析パイプラインの開発は以下のプロセスに従ってください：

### 1.1.1 要件定義フェーズ

- ユーザー（研究者・解析依頼者）との対話を通じて解析要件を明確に定義する
- 詳細な仕様書を作成し、`doc/`ディレクトリに保存する
- 仕様書には以下の要素を含める：
  - 解析の目的と生物学的背景
  - 入力データの形式と前処理状況
  - 使用するアルゴリズムとパラメータ
  - 期待される出力と結果の解釈方法
  - リソース要件（メモリ、CPU、実行時間）
  - バリデーション基準

### 1.1.2 テスト駆動開発（TDD）フェーズ

1. **テスト設計**:
   - 仕様書に基づいてテストケースを作成
   - 単体テスト、統合テスト、機能テストのレベルで計画
   - `tests/`ディレクトリに各種テストを配置

2. **テスト実装**:
   - 各機能のテストをまず実装（Red）
   - テストが意図通り失敗することを確認
   - テスト用の小規模データセットを`tests/fixtures/`に配置

3. **コード実装**:
   - テストをパスする最小限のコードを実装（Green）
   - SHIROKANEの小規模ジョブとしてテストを実行

4. **リファクタリング**:
   - コードの品質と効率を改善（Refactor）
   - パフォーマンス最適化
   - リファクタリング後もテストがパスすることを確認

5. **反復**:
   - 以上のサイクルを繰り返し、機能を段階的に拡張

### 1.1.3 仕様変更対応

- 要件変更が発生した場合は、まず仕様書を更新
- 変更に対応する新しいテストを追加
- 既存機能への影響を評価（回帰テスト実行）
- コード修正と再テスト

## 2. SHIROKANEでの環境セットアップ

### 2.1 基本的な環境設定

- moduleコマンドを使用して必要な言語環境をロードする
```bash
module load R/4.1.0    # Rの指定バージョンをロード
module load python/3.9  # Pythonの指定バージョンをロード
```

- ジョブスクリプト内では、必ず使用する言語のバージョンを明示する
- 異なる言語間で連携する場合は、互換性のあるバージョンを選択する

### 2.2 パッケージ・モジュール管理

#### R言語
- プロジェクト内のパッケージ管理には`renv`を使用する
- 新しいパッケージを追加した場合は必ず`renv::snapshot()`を実行する
- プロジェクトの初期セットアップ時に`renv::restore()`を実行する

#### Python
- 環境管理には`conda`または`venv`を使用する
- 依存関係は`requirements.txt`または`environment.yml`で明示的に管理
- 仮想環境の作成と使用：
```bash
# conda環境の作成
conda create -n myenv python=3.9 pandas numpy scipy
# 環境のアクティベート
source activate myenv
```

### 2.3 ディレクトリ構造

```
project/
├── README.md        # プロジェクト概要
├── doc/             # 要件定義書、仕様書
├── data/            # 入力データ（Git管理外）
├── results/         # 解析結果（Git管理外）
├── src/             # ソースコード
│   ├── R/           # R関数とスクリプト
│   ├── python/      # Pythonモジュールとスクリプト
│   └── bash/        # シェルスクリプト
├── job_scripts/     # ジョブスクリプト
├── notebooks/       # JupyterノートブックやRMarkdown
├── tests/           # テストスクリプト
│   ├── unit/        # ユニットテスト
│   ├── integration/ # 統合テスト
│   └── fixtures/    # テスト用データ
├── env/             # 環境設定ファイル
│   ├── renv/        # renv関連ファイル
│   ├── renv.lock    # Rパッケージバージョン情報
│   └── environment.yml # conda環境定義
└── docs/            # ドキュメント
```

## 3. ジョブスクリプトの作成と実行

### 3.1 基本的なジョブスクリプト

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -pe def_slot 4

# モジュールのロード
module load R/4.1.0
module load python/3.9

# 実行環境の情報出力
echo "Job started at $(date)"
echo "Running on host: $(hostname)"
echo "CPU cores: $NSLOTS"

# スクリプトの実行
Rscript src/R/analysis.R
python src/python/process_results.py

echo "Job finished at $(date)"
```

### 3.2 パラレル処理

#### R言語での並列処理
```r
# 並列処理のためのパッケージ
library(parallel)
library(foreach)
library(doParallel)

# 環境変数からコア数を取得
cores <- as.integer(Sys.getenv("NSLOTS", unset = "1"))

# 並列処理の設定
cl <- makeCluster(cores)
registerDoParallel(cl)

# 並列処理の実行
results <- foreach(i = 1:100, .combine = rbind) %dopar% {
  analysis_function(i)
}

# クラスターのクリーンアップ
stopCluster(cl)
```

#### Pythonでの並列処理
```python
import multiprocessing as mp
import os

# 環境変数からコア数を取得
cores = int(os.environ.get('NSLOTS', 1))

# 並列処理の実行
with mp.Pool(processes=cores) as pool:
    results = pool.map(analysis_function, range(100))
```

### 3.3 メモリ効率の高いコーディング

- 大規模データ処理では、ストリーミング処理やチャンク処理を活用する
- Rでは`data.table`や`arrow`、Pythonでは`dask`や`vaex`などのパッケージを使用
- メモリ使用量を監視し、必要に応じて中間結果をディスクに保存
- 必要に応じてガベージコレクションを明示的に呼び出す（R:`gc()`, Python:`gc.collect()`）

### 3.4 複合言語ワークフロー

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=16G
#$ -pe def_slot 8

module load R/4.1.0
module load python/3.9

# ステップ1: Bashでデータの前処理
echo "Step 1: Preprocessing data"
bash src/bash/preprocess.sh data/raw_data.txt > data/processed_data.txt

# ステップ2: Pythonで特徴抽出
echo "Step 2: Feature extraction"
python src/python/extract_features.py --input data/processed_data.txt --output data/features.csv

# ステップ3: Rで統計解析と可視化
echo "Step 3: Statistical analysis"
Rscript src/R/analyze.R --input data/features.csv --output results/

echo "Workflow completed at $(date)"
```

## 4. 言語別コーディング規約とテスト

### 4.1 R言語

- tidyverseスタイルガイドに従ってコードを記述
- 関数は一つの明確な目的を持つように設計
- roxygen2を使用したドキュメント化
- testthatによるユニットテスト実装

### 4.2 Python

- PEP 8スタイルガイドに従う
- Docstringsでの関数・クラスドキュメント化
- Type hintsを使用した型アノテーション
- pytestによるユニットテスト実装

### 4.3 Bash

- シバン（`#!/bin/bash`）を必ず記述
- エラーハンドリングを明示的に実装（`set -e`, `trap`）
- 変数はクォートで囲む
- 複雑な処理はBashではなくPythonやRに委譲

### 4.4 テスト駆動開発の実践

各言語でのTDD実践方法と具体的な例を以下に示します。

### 4.4.1 言語別テストフレームワーク

| 言語 | テストフレームワーク | ファイル配置 | 実行方法 |
|------|---------------------|-------------|---------|
| R    | testthat            | tests/testthat/ | `testthat::test_dir("tests/testthat")` |
| Python | pytest           | tests/         | `python -m pytest tests/` |
| Bash   | bats, shunit2    | tests/         | `bats tests/` |

### 4.4.2 Rのテスト例

```r
# tests/testthat/test_analysis.R
test_that("differential expression analysis works", {
  # テストデータの準備
  test_data <- load_test_data()
  test_metadata <- load_test_metadata()
  
  # 関数の実行
  results <- run_differential_expression(test_data, test_metadata, 
                                        design = ~ condition)
  
  # 結果の検証
  expect_true("gene_id" %in% colnames(results))
  expect_true("log2FoldChange" %in% colnames(results))
  expect_true("pvalue" %in% colnames(results))
  expect_true(nrow(results) > 0)
  expect_true(all(!is.na(results$pvalue)))
  
  # エッジケースのテスト
  expect_error(run_differential_expression(test_data, test_metadata[1:2,]))
})

# テスト用のジョブスクリプト
# tests/job_scripts/test_deseq2.sh
```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=4G
#$ -l mem_req=4G
#$ -pe def_slot 1

module load R/4.1.0
Rscript -e "testthat::test_file('tests/testthat/test_analysis.R')"
```

### 4.4.3 Pythonのテスト例

```python
# tests/test_analysis.py
import pytest
import pandas as pd
import numpy as np
from src.python.feature_extraction import extract_features

def test_feature_extraction():
    # テストデータの準備
    test_data = pd.DataFrame({
        'gene_id': ['gene1', 'gene2', 'gene3'],
        'sample1': [10, 20, 30],
        'sample2': [15, 25, 35]
    })
    
    # 関数の実行
    features = extract_features(test_data)
    
    # 結果の検証
    assert "gene_id" in features.columns
    assert "mean_expression" in features.columns
    assert "variance" in features.columns
    assert len(features) == 3
    
    # エッジケースのテスト
    with pytest.raises(ValueError):
        extract_features(pd.DataFrame())
```

### 4.4.4 Bashのテスト例

```bash
# tests/test_preprocessing.bats
#!/usr/bin/env bats

setup() {
    # テスト環境の準備
    export TEST_DIR="$(mktemp -d)"
    cp tests/fixtures/sample.fastq "$TEST_DIR/"
}

teardown() {
    # テスト環境のクリーンアップ
    rm -rf "$TEST_DIR"
}

@test "FastQC quality control works" {
    # スクリプトの実行
    bash src/bash/run_fastqc.sh "$TEST_DIR/sample.fastq" "$TEST_DIR/output"
    
    # 結果の検証
    [ -f "$TEST_DIR/output/sample_fastqc.html" ]
    [ -f "$TEST_DIR/output/sample_fastqc.zip" ]
}
```

### 4.4.5 統合テスト例

複数の言語やステップを組み合わせた統合テストの例:

```bash
# tests/integration/test_full_pipeline.sh
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -pe def_slot 2

set -e

echo "Running integration test for full pipeline"

# テスト環境の準備
TEST_DIR="$(mktemp -d)"
cp -r tests/fixtures/* "$TEST_DIR/"

# ステップ1: Bash前処理
echo "Step 1: Running preprocessing"
bash src/bash/preprocess.sh "$TEST_DIR/raw_data.txt" > "$TEST_DIR/processed_data.txt"

# ステップ2: Python特徴抽出
echo "Step 2: Running feature extraction"
module load python/3.9
python src/python/extract_features.py --input "$TEST_DIR/processed_data.txt" --output "$TEST_DIR/features.csv"

# ステップ3: R解析
echo "Step 3: Running statistical analysis"
module load R/4.1.0
Rscript src/R/analyze.R --input "$TEST_DIR/features.csv" --output "$TEST_DIR/results/"

# 結果の検証
echo "Verifying results"
[ -f "$TEST_DIR/results/differential_expression.csv" ]
[ -f "$TEST_DIR/results/pca_plot.pdf" ]
[ -f "$TEST_DIR/results/heatmap.pdf" ]

# テストレポート生成
echo "Generating test report"
Rscript src/R/generate_test_report.R --input "$TEST_DIR/results/" --output "$TEST_DIR/test_report.html"

echo "Integration test completed successfully"
```

## 5. デバッグとロギング

### 5.1 基本方針

- すべてのスクリプトにログ機能を実装
- エラー時には具体的なメッセージを出力
- 長時間実行ジョブには進捗状況の出力を含める

### 5.2 言語別ロギング実装

#### R言語
```r
library(futile.logger)
flog.threshold(INFO)
log_file <- file.path("logs", paste0("analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
flog.appender(appender.file(log_file))

flog.info("Analysis started")
# 処理
flog.info("Analysis completed")
```

#### Python
```python
import logging
logging.basicConfig(
    filename=f"logs/analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

logger.info("Analysis started")
# 処理
logger.info("Analysis completed")
```

#### Bash
```bash
log_file="logs/script_$(date +%Y%m%d_%H%M%S).log"
mkdir -p $(dirname $log_file)

log() {
  echo "[$(date +%Y-%m-%d\ %H:%M:%S)] $1" | tee -a $log_file
}

log "Script started"
# 処理
log "Script completed"
```

## 6. データ管理と再現性

### 6.1 データの保存と管理

- 生データは読み取り専用として扱い、処理済みデータと明確に区別
- 中間ファイルは適切に名前付けし、生成プロセスを記録
- 大規模データはSHIROKANEの適切なストレージ領域に保存

### 6.2 再現性の確保

- すべての乱数シードを明示的に設定（R:`set.seed()`, Python:`np.random.seed()`）
- 環境依存性を最小限に抑え、必要な依存関係を明示的に記録
- 解析パラメータはYAMLやJSONなどの設定ファイルに保存
- ワークフローの各ステップを再実行可能な形で設計

### 6.3 バージョン管理

- Gitを使用したコード管理
- データやモデルはGit LFS、DVC、またはS3などで管理
- コミットメッセージには変更内容を明確に記述
- ブランチ戦略を定義し、一貫して適用

## 7. 言語間連携のベストプラクティス

### 7.1 データ交換形式

- CSVやTSVなどの標準テキストフォーマットを基本とする
- 大規模データにはParquet、HDF5、SQLiteなどを使用
- JSONやYAMLは設定や小～中規模データに使用

### 7.2 プロセス間通信

- パイプとリダイレクトを活用した基本的なデータフロー
```bash
python generate_data.py | R --vanilla -f process_data.R > results.txt
```

- REPLベースの連携（R/Python）
```r
# RからPythonコードを実行
library(reticulate)
py_run_string("import numpy as np; data = np.random.randn(100, 5)")
r_data <- py$data
```

### 7.3 ワークフロー管理

- 基本的なMakefileによる依存関係管理
- Snakemakeなどのワークフローマネージャの活用
- ステップごとの入出力を明確に定義

## 8. SHIROKANEリソース最適化

### 8.1 リソース要求の最適化

- 必要に応じてジョブクラスを選択（共有メモリ、GPUなど）
- メモリ要求を適切に設定（`s_vmem`, `mem_req`）
- 並列コア数を処理内容に応じて調整（`pe def_slot`）

### 8.2 ストレージ利用の最適化

- ストレージ階層を意識したデータ配置（ホーム、グループ、一時領域）
- 大規模一時ファイルには`/work`領域を使用
- 長期保存データには適切なバックアップ戦略を実装

### 8.3 ジョブチェーン実行

```bash
# ステップ1: データ前処理
JOB1=$(qsub job_scripts/01_preprocess.sh)

# ステップ2: 主要解析（前処理完了後に実行）
JOB2=$(qsub -hold_jid $JOB1 job_scripts/02_analysis.sh)

# ステップ3: 結果の可視化（主要解析完了後に実行）
JOB3=$(qsub -hold_jid $JOB2 job_scripts/03_visualization.sh)
```

## 9. 言語選択のガイドライン

特定のタスクに最適な言語を選択するためのガイドライン：

### 9.1 R言語が適しているケース

- 統計解析と可視化が主要なタスク
- Bioconductorパッケージを活用した生命科学解析
- インタラクティブな探索的データ分析
- データフレーム操作や機械学習モデル評価

### 9.2 Pythonが適しているケース

- カスタム機械学習モデルの開発と訓練
- データ前処理パイプラインの構築
- 大規模データの並列処理
- 複雑なデータ構造の操作やWeb API連携

### 9.3 Bashが適しているケース

- シンプルなファイル操作やデータ前処理
- 既存コマンドラインツールの連携
- ジョブスケジューリングと監視
- システムレベルのリソース管理

### 9.4 複合アプローチ

- 各言語の強みを生かした複合ワークフロー
- 処理速度とメモリ効率のバランスを考慮
- 開発速度とコード保守性のトレードオフを評価

## 10. 生物情報学固有のルール

- バイオインフォマティクスツールの適切なバージョン選択
- 参照ゲノムとアノテーションの明示的な指定
- 解析パラメータの妥当性確認と文書化
- 生物学的解釈に基づく結果の検証と評価 