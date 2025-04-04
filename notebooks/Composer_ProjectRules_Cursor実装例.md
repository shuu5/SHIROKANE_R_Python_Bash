# SHIROKANE R/Python/Bash のための Cursor ルールファイル

#cursor #composer #rules_for_ai #R #Python #Bash #スパコン #SHIROKANE #ジョブスクリプト #AGE #バージョン管理

> [!NOTE]
> SHIROKANEスーパーコンピュータでのR/Python/Bash言語解析に関するCursorのルールファイルです。

## Cursorでのプロジェクトルール設定方法

Cursorでプロジェクトルールを適用するには、`.cursor/rules`ディレクトリに`.mdc`拡張子のルールファイルを作成します。

### ルールファイルの設置

1. プロジェクトのルートディレクトリに`.cursor/rules`ディレクトリを作成
2. 任意の名前で`.mdc`拡張子のファイルを作成

### ルールファイルの内容

```
---
Description: SHIROKANE スーパーコンピュータでのR/Python/Bashによる生物情報学解析ガイド
Globs: **/*.R, **/*.py, **/*.sh, **/*.Rmd, **/*.ipynb, job_scripts/**/*
---

# SHIROKANE R/Python/Bash 解析のための Project Rules

## 知識ベース

- SHIROKANE は東京大学ヒューマンゲノムセンターが提供するスーパーコンピュータシステム
- 計算ノードで実行するジョブはジョブスケジューラ（UGE）を通じて投入
- moduleコマンドで言語環境やツールをロード（例: `module load R/4.1.0`）
- 並列処理のために SGE の環境変数 `$NSLOTS` を利用

## 基本原則

1. **再現性を最優先** - 処理手順の明示的スクリプト化とバージョン情報の記録
2. **適切なリソース管理** - ジョブスクリプトでの適切なリソース要求
3. **言語間連携** - R、Python、Bashの強みを活かした設計
4. **エラーハンドリング** - 堅牢なエラー処理と詳細なログ出力
5. **データ管理** - 適切なデータ保存とバックアップ戦略
6. **要件駆動開発** - 要件定義と仕様書作成を経た実装
7. **テスト駆動開発** - テストを先に書き、要件を満たすコード実装

## 要件定義と開発プロセス

1. **要件定義**: ユーザーとの対話による機能や制約条件の明確化、詳細な仕様書作成
2. **テスト駆動開発（TDD）**: 仕様書に基づいたテストケース作成後に実装
3. **仕様変更の管理**: 要件変更発生時は仕様書更新→新テスト作成→コード修正の順

## ジョブスクリプトの基本構造

- 適切なシェバン（`#!/bin/bash`）
- AGEディレクティブ（`#$`で始まる行）
- カレントディレクトリ指定（`#$ -cwd`）
- リソース指定（メモリ、CPUコア数）
- モジュールロード
- 実行環境情報の出力
- 処理開始・終了のログ出力

## 言語ごとの最適な使用ケース

- **R**: 統計解析、可視化、Bioconductorパッケージを用いた解析
- **Python**: 機械学習、大規模データ処理、カスタムパイプライン
- **Bash**: ファイル操作、ワークフロー管理、既存ツールの連携

## 言語間のデータ連携

1. 標準的なファイルフォーマット（CSV、TSV）による受け渡し
2. 大規模データには効率的なフォーマット（Parquet、HDF5）を使用
3. R-Python連携は reticulate または rpy2 を活用

## パラレル処理の実装

**R言語**: `parallel`パッケージ使用、`NSLOTS`からコア数取得、`stopCluster()`で解放
**Python**: `multiprocessing`モジュール使用、Context manager活用
**Bash**: GNU Parallel または xargsによる並列処理

## エラーハンドリングとロギング

**R**: `futile.logger`パッケージ、日時情報含むログ、適切なエラー終了
**Python**: 標準`logging`モジュール、try-except捕捉
**Bash**: `set -e`と`trap`活用、ログ関数定義

## テスト実装

**ユニットテスト**: R(`testthat`)、Python(`unittest`)
**統合テスト**: ワークフロー全体のテスト、テストデータセット準備

## メモリ効率の最適化

1. チャンク処理とストリーミング処理
2. メモリ使用量の監視と管理
3. 適切なパッケージ選択(R: data.table/arrow, Python: dask/vaex)

## ジョブチェーン実行

依存関係のあるジョブ実行には`-hold_jid`オプション使用、ワークフロー制御スクリプト準備

## 推奨するディレクトリ構造

- `doc/`: 仕様書や設計ドキュメント
- `data/`, `results/`: 入力データと解析結果（Git管理外）
- `src/`: ソースコード（言語別にサブディレクトリ分け）
- `job_scripts/`: ジョブスクリプト
- `tests/`: テストスクリプト（unit/integration/test_data）
- `notebooks/`, `env/`, `logs/`: ノートブック、環境設定、ログ

## コーディング規約

- R: tidyverseスタイルガイド
- Python: PEP 8
- Bash: Google Shell Style Guide

## 参考リンク

- SHIROKANE ユーザーガイド: https://gc.hgc.jp/uge/
- R言語公式ドキュメント: https://www.r-project.org/
- Python公式ドキュメント: https://www.python.org/doc/
- Bioconductorパッケージ: https://www.bioconductor.org/
- テストフレームワーク: 
  - R: https://testthat.r-lib.org/
  - Python: https://docs.python.org/3/library/unittest.html
- TDDの原則: https://www.agilealliance.org/glossary/tdd/ 