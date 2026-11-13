CDC clinical：Nanopore 單樣本病原體探索

資料與樣本
  samples.tsv：四個來源檔案的追蹤表；sample_001 為暫定生物樣本識別碼。
  barcode17 / NB17 / SG26-GS01：DNA library。
  barcode18 / NB18 / ST26-UM01：cDNA library。
  兩種 library 以 equal-fmole 混合，同一份 pool 分兩次定序，屬於技術重複。
  檔名前綴 PBE34168、PBI71026 用於追蹤兩批資料；尚未假定與 PS010/PS011 的對應。
  檢體類型與陰性對照資訊尚未提供，詳見 sample_metadata.json。

環境建置
  bash scripts/setup_environment.sh
  獨立環境：.conda/envs/nanopore/；不變更 base 環境。
  environment.yml 為相依工具規格；environment-linux-64.lock.txt 為完整下載版本鎖定。
  需要重建完全相同版本時，使用：
    conda create -p .conda/envs/nanopore --file environment-linux-64.lock.txt

來源掃描
  python scripts/nanopore_pipeline.py scan --manifest samples-new.tsv
  遞迴掃描 fastq_pass/，遇到未知 barcode、空目錄或重複來源會停止。
  已有的樣本表不會被掃描覆寫；執行前會驗證檔案存在與檔案大小。

小份驗證
  bash batch_dehost.sh --limit-reads 1000 --output validation/smoke --keep-qc
  每個來源取前 1000 reads，總計 4000 reads；並非隨機抽樣，不能代表病原體組成。
  已知答案的控制測試：
    PATH="$PWD/.conda/envs/nanopore/bin:$PATH" python scripts/validate_controls.py
  控制組包含人類 genomic／合成剪接 reads、lambda 噬菌體、低品質與過短 reads。
  lambda 只用於流程測試；不代表檢體中有此生物。

全量去宿主（需另外啟動）
  bash batch_dehost.sh --output clean_results --threads 16
  預設 Q>=10、最短長度 200 bp、頭尾各裁 15 bp，沿用初版設定。
  固定裁切不等於完整 adapter / primer removal；需根據建庫與 basecalling 紀錄確認。
  可設 --headcrop 0 --tailcrop 0 關閉固定裁切，或調整 --quality / --min-length。
  原始 FASTQ 不會改動。每批先 QC、去宿主，再按 library 合併技術重複。
  DNA 使用 map-ont，cDNA 使用 splice preset，以處理跨剪接位點的人類轉錄序列。
  使用官方 human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401 參考。
  下載後以官方 SHA256 核對；DNA/cDNA 各建一份 minimap2 索引。
  此為 masked 參考，目的是減少微生物誤刪；nonhost 標籤不保證人類序列完全移除。

結果
  sources/<source_id>/：逐批 clean.fastq.gz、chopper.log、hostile.log、hostile.json。
  source_summary.tsv：逐批 raw / QC / nonhost reads、bases 與去宿主數量。
  sample_summary.json：參數、來源與分 library 合併結果。
  sample_001.SG26-GS01.nonhost.fastq.gz：DNA 的兩次定序合併。
  sample_001.ST26-UM01.nonhost.fastq.gz：cDNA 的兩次定序合併。
  合併後重新統計，確認 reads 與 bases 等於各來源總和。
  預設刪除 QC 暫存，使用 --keep-qc 保留。
  同一輸出目錄有檔案鎖；同參數可續跑。來源或參數更動時，請使用新輸出目錄。
  合併 gzip 使用 concatenated gzip members；seqkit/pigz 支援此標準格式。

後續病原體探索
  已安裝 Kraken2 與 BLAST；尚未下載廣泛分類資料庫，也尚未執行全量病原體分類。
  先選用有日期及 checksum 的資料庫，涵蓋細菌、病毒、真菌與原蟲；病毒庫或細菌庫
  單獨使用都不足以回答開放式的病原體問題。
  DNA 與 cDNA 分別分類並保留兩次上機的批次結果，再整合候選物種。
  候選物種須回貼參考基因組，評估獨立 reads、分布、覆蓋廣度／深度、比對品質，
  並與人類及近緣物種做競爭比對，避免把局部保守序列判成物種感染。
  equal-fmole 不等於相同 read 數或相同微生物量，不直接比較 DNA/cDNA 的分類百分比。
  沒有陰性對照時，低量背景、試劑污染、barcode 串擾與真實訊號仍需分辨。
  分類命中為候選微生物證據；是否為致病原因需結合檢體來源、臨床背景與獨立確認。

官方文件
  https://github.com/bede/hostile
  https://github.com/wdecoster/chopper
  https://github.com/lh3/minimap2
  https://github.com/DerrickWood/kraken2/wiki/Manual

2026-11-12 全量執行
  去宿主：bash batch_dehost.sh --output clean_results --threads 38 --jobs 4
    --jobs 為新增選項，4 個來源平行處理（大檔先跑）；與逐一處理的小份結果逐欄相同。
    log：logs/full-dehost.log
  Kraken2 資料庫：Kraken2 PlusPF 2026-10-22（archaea、bacteria、viral、plasmid、human、
    UniVec_Core、protozoa、fungi），來源 https://benlangmead.github.io/aws-indexes/k2
    databases/kraken2/k2_pluspf_20261022/；以檔案大小與 S3 multipart ETag 核對（官方未提供 md5）。
    scripts/prepare_kraken2_db.sh（核對＋解壓），log：logs/kraken2-db-download.log
  分類：scripts/run_kraken2.sh（逐批來源與合併 library 各一份 report，含 minimizer 資料）
    彙整：clean_results/kraken2/overview.tsv、species_table.tsv、genus_table.tsv
  自動串接：scripts/auto_kraken2_after_dehost.sh，log：logs/auto-kraken2.log
