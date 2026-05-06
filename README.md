# bam_stat

BAM 对齐一致性（Concordance）统计工具。

## 简介

对 BAM 文件中每条 Read 的 CIGAR 对齐信息逐碱基分析，排除参考区域（HC Regions）和已知变异位点（HC Variants）的干扰，输出每条 Read 的匹配/错配/插入/删除等统计指标及一致性评分。

## 用法

```bash
# 构建
cargo build

# 运行
cargo run -- bam-concordance <reffasta> <aligned_bam> [OPTIONS]

# 参数
  --hcregions <BED>    一致性计算排除的区域（BED 格式）
  --hcvariants <VCF>   一致性计算排除的已知变异位点（VCF 格式）
  --chrom <CHROM>      染色体过滤（当前未实现）
```

输出文件：`{aligned_bam 文件名前缀}.metric.csv`

## 统计逻辑

### 整体流程

程序使用**生产者-消费者多线程**管道处理 BAM 文件：

1. **预加载**（主线程）：在独立线程中分别解析 BED 区域文件和 VCF 变异文件，再将完整参考 FASTA 加载进内存
2. **BAM 读取**：10 个线程读取 BAM，每条记录提取关键字段包装成 `RecordReplica`，通过 `crossbeam::channel`（容量 200）发送给消费者
3. **统计计算**：8 个工作线程从 channel 拉取记录，逐 CIGAR 操作与参考序列比对，产出 `Stat`
4. **结果输出**：主线程接收所有 Stat，以 TSV 格式写入 CSV 文件

### 逐 CIGAR 操作统计

`stat_record_core()` 是核心函数，维护两个游标：
- **`rpos_cursor`** — 参考序列上的当前位置
- **`qpos_cursor`** — 查询序列（Read）上的当前位置

遍历 CIGAR 的每个操作，按以下规则统计：

#### 1. 匹配与错配（M / MM）

| CIGAR 操作 | 统计字段 | 含义 |
|-----------|---------|------|
| `Equal(N)` | `stat.m += N` | 精确匹配的碱基数 |
| `Diff(N)` | `stat.mm += N` | 替换错配的碱基数 |

**有效性检查**：只有当该位置同时满足 "在 HC 区域内" 且 "不在 VCF 变异位点" 时，才计入 m/mm。排除的碱基数累加到 `ignore_bps`。同时 `Diff` 会额外记录错配的参考位置列表 `mm_ref_positions`，用于后续分析。

#### 2. 插入（Insertion）

| CIGAR 操作 | 统计字段 | 含义 |
|-----------|---------|------|
| `Ins(N)` | `stat.h_ins` | 同源多聚体（Homopolymer）插入碱基数 |
| `Ins(N)` | `stat.non_h_ins` | 非同源多聚体插入碱基数 |

HP 判定：比较插入的碱基序列与参考序列在插入位置前一个碱基和后一个碱基，如果完全一致则判定为同源多聚体插入。

```
valid_insertion(ref_pos, ref_name) → 检查 ref_pos 及 ref_pos-1 是否在 HC 区域且非变异位点
    └─ 如果在：
        检查插入序列 cur_bp_seq 是否与邻近参考碱基相同
        相同 → stat.h_ins += N    （HP 插入）
        不同 → stat.non_h_ins += N（非 HP 插入）
        同时记录 stat.ins_ref_positions
    └─ 如果不在：
        stat.ignore_bps += N
```

#### 3. 删除（Deletion）

| CIGAR 操作 | 统计字段 | 含义 |
|-----------|---------|------|
| `Del(N)` | `stat.h_del` | 同源多聚体（Homopolymer）删除碱基数 |
| `Del(N)` | `stat.non_h_del` | 非同源多聚体删除碱基数 |

HP 判定：对被删除的参考碱基，检查其左右邻居碱基是否相同（同源多聚体特征）。

```
valid_point(ref_pos, ref_name) → 检查是否在 HC 区域且非变异位点
    └─ 如果在：
        检查该位置参考碱基与其左右邻居是否相同
        相同 → stat.h_del += 1     （HP 删除）
        不同 → stat.non_h_del += 1（非 HP 删除）
        同时记录 stat.del_ref_positions
    └─ 如果不在：
        stat.ignore_bps += 1
```

#### 4. 其他 CIGAR 操作

`SoftClip`、`HardClip`、`RefSkip(N)` 等不产生统计贡献，仅推进游标。

游标推进规则：

| CIGAR 操作 | rpos_cursor 变化 | qpos_cursor 变化 |
|-----------|-----------------|-----------------|
| `Equal/Diff(N)` | += N | += N |
| `Del/RefSkip(N)` | += N | 无变化 |
| `Ins/SoftClip(N)` | 无变化 | += N |

### 统计指标

Stat 结构体的核心字段：

| 字段 | 含义 |
|------|------|
| `ch` | 通道/频道编号（来自 BAM aux 标签 `ch`） |
| `q_len` | Read 序列长度 |
| `passes` | subread passes（来自 BAM aux 标签 `np`） |
| `m` | 匹配碱基数（Equal） |
| `mm` | 错配碱基数（Diff） |
| `h_ins` | 同源多聚体插入碱基数 |
| `non_h_ins` | 非同源多聚体插入碱基数 |
| `h_del` | 同源多聚体删除碱基数 |
| `non_h_del` | 非同源多聚体删除碱基数 |
| `ignore_bps` | 被排除的碱基数 |
| `ref_start` / `ref_end` | 比对在参考序列上的起止位置 |
| `rq` | Read quality（来自 BAM aux 标签 `rq`） |

### 派生指标

| 指标 | 计算方式 | 含义 |
|------|---------|------|
| `ins_bp()` | `non_h_ins + h_ins` | 总插入碱基数 |
| `del_bp()` | `non_h_del + h_del` | 总删除碱基数 |
| `align_span()` | `ins_bp + del_bp + m + mm` | 对齐跨度（含 indel 的总覆盖碱基数） |
| `concordance()` | `m / align_span` | **一致性率**：匹配碱基占对齐跨度的比例 |
| `concordance_qv()` | `-10 * log10(1 - concordance)` | 一致性 Phred 质量值，错误率下限 1e-6（对应 Q60） |
| `query_converage()` | `(ins_bp + m + mm) / q_len` | 查询序列覆盖率：被计入对齐的查询碱基数 / 序列总长 |
| `predictedConcordance` | `rq` 字段的值 | Read 自身报告的预测一致性（来自 aux 标签 `rq`） |

### HC（Hide / Exclude from Concordance）机制

两条独立的过滤规则：

1. **HC Regions（BED）**：用户通过 `--hcregions` 指定的区域，这些区域内的碱基不计入一致性统计（用于排除低质量比对区域）
2. **HC Variants（VCF）**：用户通过 `--hcvariants` 指定的已知变异位点，这些位点及其对应的同源多聚体相邻位点被排除，避免已知杂合/同源多聚体位点污染一致性评分

两个条件用 `&=` 连接——碱基必须同时满足"在 HC 区域内"且"不在变异位点"才算有效。

### 输出字段

TSV 格式，每行一条 Read：

| 列名 | 含义 |
|------|------|
| `channel_id` | 通道编号 |
| `readLengthBp` | Read 长度 |
| `subreadPasses` | subread passes |
| `queryConverage` | 查询覆盖率 |
| `predictedConcordance` | 预测一致性（rq 标签） |
| `concordance` | 一致性率 |
| `concordanceQv` | 一致性 Phred 质量值 |
| `matchBp` | 匹配碱基数 |
| `mismatchBp` | 错配碱基数 |
| `nonHpInsertionBp` | 非 HP 插入碱基数 |
| `nonHpDeletionBp` | 非 HP 删除碱基数 |
| `hpInsertionBp` | HP 插入碱基数 |
| `hpDeletionBp` | HP 删除碱基数 |
| `ignoreBp` | 被排除的碱基数 |
| `refStart` | 参考起始位置 |
| `refEnd` | 参考结束位置 |
| `mmRefPositions` | 错配参考位置（逗号分隔字符串） |
| `insRefPositions` | 插入参考位置（逗号分隔字符串） |
| `delRefPositions` | 删除参考位置（逗号分隔字符串） |
