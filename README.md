# probebwa

[![CI](https://github.com/GUIBA-EX/probebwa/actions/workflows/ci.yml/badge.svg)](https://github.com/GUIBA-EX/probebwa/actions/workflows/ci.yml)

UCE(超保守元件)探针/bait 序列比对到参考基因组的短读比对器:密集 k-mer 种子
扩展、banded affine-gap 比对(质量感知错配罚分 + 同聚物感知 gap-open)、贝叶斯
后验 MAPQ 打分;支持单端和配对末端(候选 shortlist + 联合后验,能把落进重复
区域的读段纠正到正确拷贝)。

当前正在按 **Stampy 1.0.32 (r3761) 输出优先、之后再加速** 的标准修正。
尚未达到完整输出一致，不能将现版本视为 Stampy 的等价替代品。
原版计算结果覆盖随机匹配统计、漏检概率及片段长度更新；候选搜索、比对
计分/二次精细比对、比对路径置信度、配对细节和 SAM 标签仍需逐项还原。

## 编译

```bash
cargo build --release
```

默认按本机 CPU 特性编译(`-C target-cpu=native`);要编译一次、拷贝到别的机器
上跑,用更保守的基线覆盖:

```bash
RUSTFLAGS="-C target-cpu=x86-64-v2" cargo build --release
```

## 用法

```bash
# 建索引
./target/release/probebwa build-genome --species human --assembly hg38 -G hg38 ref.fa.gz
./target/release/probebwa build-hash --genome hg38 -H hg38

# 比对(两个 read 文件即视为配对末端)
./target/release/probebwa map --genome hg38 --hash hg38 -M reads_1.fq.gz reads_2.fq.gz > out.sam

# BAM 输出 + 坐标排序索引
./target/release/probebwa map --genome hg38 --hash hg38 \
    --outputformat bam --output out.bam --index -M probes.fq.gz
```

完整参数(gap 罚分、结构变异先验、多线程、接头剪切、CASAVA 1.8 兼容、
read group 等)见 `probebwa map --help`。

初始片段长度可通过 `--insertsize 250 --insertsd 80` 指定。SD 默认值 80
与 Stampy 1.0.32 可执行实现一致（其旧帮助文字写的是 60）；与旧版
probebwa 或显式使用 SD=60 的实验比较时，请明确加上 `--insertsd 60`。

## 架构

| 模块 | 作用 |
|---|---|
| `index/` | 基因组索引(`.stidx`)+ 15-mer 哈希表(`.sthash`) |
| `mapper/` | 种子扫描、候选聚类、单端/配对末端比对 |
| `align/` | Banded affine-gap 比对 |
| `mapq/` | 贝叶斯后验 MAPQ、插入片段大小在线学习 |
| `io/` | FASTA/FASTQ 解析、SAM/BAM 输出 |

## 测试

```bash
cargo test --release
```

`.github/workflows/ci.yml` 在每次 push/PR 到 `main` 时自动跑 build、test、
`clippy --all-targets -- -D warnings`。

原版差分检查（任一比对字段、标签、记录数量/顺序或非运行元数据头部
不一致即返回非零；不使用坐标容差）：

```bash
python3 scripts/check_stampy_fidelity.py stampy.sam probebwa.sam --report fidelity.json
python3 -m unittest discover -s scripts -p 'test_*.py'
```

检查仅统一 `RNEXT` 的 `=` 写法和标签排列顺序，`@PG`/`@CO` 运行元数据
另列于报告。统计组件的三个 TSV 文件保存 1,480 个原版计算结果，
普通 Rust 测试无需安装 Stampy；重新生成时需要本地原版及其 Python 2.7：

```bash
/path/to/python2.7 scripts/stampy_stats_oracle.py /path/to/stampy tests/fixtures
```

这些组件测试通过不代表完整比对器一致。发布“Stampy 等价”版本前，
完整 SAM 差分检查也必须通过；现阶段该检查仍会报告差异。

### 分阶段还原进度

生产流程已接入 `align/sensitive.rs`、`align/stampy_band.rs` 与
`align/confidence.rs`。候选先做粗比对，再按相对得分筛选、去重及精细
比对；配对输出使用 BAQ 推导的比对置信度。组件通过率仍不等于最终
SAM 一致率，尤其是重复区选点、救援搜索和配对慢路径。

| 部分 | 当前验证 |
|---|---|
| 原版调用追踪 | 候选、精细比对、配对输出；追踪前后 2,000 条 SAM 记录一致 |
| 窗口选择所用粗比对 | 600 个原版样例，得分、路径、起点一致 |
| 无带宽精细比对 | 305 个合成样例；另回放 2,630 次真实调用，一致 |
| 固定带宽精细比对 | 305 个原版样例，一致 |
| 原版窗口选择 | 2,630 次真实调用一致 |
| 窗口选择与精细比对串联 | 2,630 次真实调用的得分、起点、主/备路径一致，含阈值退出 |
| BAQ/比对置信度 | 209 次原版调用；主/备质量字符串及 BAQ 修正值逐字节一致，概率误差 ≤ 1e-10 |
| 主流程 | 已接入；单端/配对精细比对分别使用调整后/原始质量值 |
| SAM 输出 | 补充 UQ、NM、PQ、SM、MQ、XQ 及参考元数据；未比对 mate 复制坐标 |
| 尚未完成 | 重复区选点、救援及配对慢路径的完整等价、最终 SAM 零差异 |

真实调用来自本地大肠杆菌样本的前 1,000 对读段，原版为 Stampy
1.0.32，单线程，片段长度初始参数 250/60。回归文件包含调用输入，
不依赖本机原版或外部基因组；两个较大文件以确定性 gzip 保存。
这些是固定样本上的精确核对，不构成对所有输入的等价证明。

原版观察工具（需单线程，完整参数透传给原版）：

```bash
/path/to/python2.7 scripts/trace_stampy.py --trace native.jsonl /path/to/stampy.py -- \
  --threads=1 -g genome -h genome -o native.sam -M reads1.fq,reads2.fq
python3 scripts/summarize_stampy_trace.py native.jsonl --report stages.json
```

Rust 候选追踪通过 `RUST_LOG=probebwa::fidelity=trace` 启用，默认关闭。
原版追踪中的坐标保留其有符号基因组坐标，Rust 候选坐标是染色体内
零基坐标，不能直接逐值比较两者。

重新生成精细比对回放输入并检查：

```bash
/path/to/python2.7 scripts/stampy_sensitive_oracle.py /path/to/stampy unbanded.tsv native.jsonl
/path/to/python2.7 scripts/stampy_window_oracle.py /path/to/stampy native.jsonl windows.tsv
python3 scripts/extract_sensitive_trace.py native.jsonl sensitive.tsv
STAMPY_SENSITIVE_FIXTURE=unbanded.tsv STAMPY_WINDOW_FIXTURE=windows.tsv \
  STAMPY_FULL_SENSITIVE_FIXTURE=sensitive.tsv cargo test --release --test sensitive
```

串联回放核对得分与比对路径；独立的 `stampy_confidence.tsv.gz` 核对
BAQ 及后验置信度。默认 SAM 与原版一样不输出 BQ/YQ 字符串，但配对
置信度已使用对应证据。原版退化窗口中的有限零概率哨兵及 x86 整数
转换行为也被显式保留。哈希索引格式、构建功能及种子查找代码未修改。

### 全流程差分检查（40,000 条记录）

同一大肠杆菌参考、20,000 对 100 bp 读段，片段长度初始参数 250/60，
与原版 Stampy 1.0.32 (r3761) 逐字段比较。“修正前”指此次配对置信度
修正前的开发快照，不是 v0.1.2 release。

| 指标 | 修正前 | 当前 |
|---|---:|---:|
| 11 个核心字段全部一致 | 39,551 | 39,565（98.9125%） |
| 核心字段及所有标签一致 | 38,554 | 38,568（96.42%） |
| FLAG 不一致 | 62 | 62 |
| POS / PNEXT 不一致（各自） | 54 | 54 |
| MAPQ 不一致 | 426 | 412 |
| CIGAR 不一致 | 21 | 21 |
| TLEN 不一致 | 80 | 80 |
| SEQ 不一致 | 25 | 25 |
| 标签不一致 | 1,283 | 1,265 |

RNAME 和 RNEXT 各有 6 条差异。记录数量、顺序、非运行元数据头部、
QNAME 和 QUAL 均一致。1/4 线程的全部 40,000 条记录和标签一致。
字段差异会重叠，不能相加。剩余 25 条 SEQ 差异全部伴随链方向变化，
其中 3 条还涉及映射状态差异；不是通过修改 SAM 序列字段来消除差异。

修正内容：反向链按带方向的坐标顺序聚类，粗比对将原始读段对齐到
相应方向的参考序列；保留粗比对坐标和路径参与候选排序，避免精比对
移动坐标后改变同分优先级。配对同分按原版候选元组的分数、读段路径、
参考路径及带方向的坐标比较；最大似然选择保留首选候选。整数配对分数
先相加再转为概率，避免浮点舍入人为打破同分。独立回退也使用粗比对
分数，与正常选点保持一致。

配对置信度进一步对齐：是否偏离单端首选位置，比较原始粗比对锚点
与最终带方向的坐标；配对选择和 PQ 使用精比对前的候选坐标，输出
片段长度与在线统计仍使用最终坐标。相同带方向坐标上的重复读段按
精比对分数保留较优者，同分保留第二端，并回退为单端置信度。
反向读段按 5′ 端判定重合，不能只比较 SAM 左端位置。


**完整 SAM 严格差分仍失败，不能发布为 100% Stampy 等价实现。**
本轮 PQ 差异 924 → 920，SM 和 MQ 各 497 → 488。候选首选与
原版仍不完全相同，会影响配对置信度及在线片段长度统计，需继续核对。
原版同分轮换计数器也尚未复现：直接观测 40,000 次候选调用，原版有
37,414 次推进 2 步、2,524 次推进 3 步、62 次不推进。用当前候选结果
推算有 67 次推进差异，试验令 SEQ 差异从 25 回升到 37，故未保留。
不能在候选枚举和处理阶段尚未对齐时，用猜测调用次数替代原版状态。

当前 x86-64 release 测试 55 项、Python 检查器测试 6 项及 Clippy 通过。
新增前 1,000 对读段中 16 组、93 个原版同分配对候选的排序回放；反向重复区端到端用例
已另行运行原版确认，选择读段 1 的反向位置 901 和读段 2 的正向位置
601。反向聚类边界和最大似然同分稳定性也有回归测试。新增两个端到端
回归：正反向重复读段退回单端报告（原版存活端 MAPQ=96），以及精
比对移动锚点后将单端 SM 判为 0。后一个缩小参考用例只核对位置、
CIGAR 和 SM/MQ，不代表该用例的完整配对后验已复现。

排序样例可重建：

```bash
python2 scripts/trace_stampy.py --trace native.jsonl --candidates-only \
  --pair-candidates /path/to/stampy.py -- [Stampy mapping arguments]
python3 scripts/extract_pair_ties.py native.jsonl stampy_pair_ties.tsv.gz
```

`--tie-counter` 可在保留 ELF 符号的原版扩展上只读观测轮换计数器；
此观察器的完整 SAM 已与不加观察器的原版核对一致。使用单个原版工作
进程采样；计数器和日志顺序依赖该进程。哈希代码、索引格式及种子查找
代码未修改。剩余工作集中于候选枚举、轮换状态、配对救援及后验标签。

### 性能优化（2026-10-11）

x86-64 上的粗比对改为 SSE2 八通道计算，其他架构保留标量实现。
精比对矩阵按参考列紧凑保存实际访问的搜索带，带外仍返回原版有限
哨兵值；计算顺序、整数回绕、同分规则和主/备回溯保持不变。

AMD Ryzen 9 9950X、Linux x86-64、release 构建；使用上述 20,000 对
100 bp 读段和 250/60 初始参数。优化前后交替运行各三次，取中位数；
计时包含索引读取和 SAM 写出。这里的“优化前”是本次性能改动前的
开发构建，不是 v0.1.2，也不是原版 Stampy。

| 线程 | 优化前耗时 | 优化后耗时 | 耗时减少 | 优化前峰值 RSS | 优化后峰值 RSS |
|---|---:|---:|---:|---:|---:|
| 1 | 12.50 秒 | 7.34 秒 | 41.3% | 141.0 MiB | 129.7 MiB |
| 4 | 5.08 秒 | 2.87 秒 | 43.5% | 138.3 MiB | 142.8 MiB |

单线程三轮：优化前 12.46/12.50/12.56 秒，优化后 7.34/7.34/7.42 秒。
四线程三轮：优化前 5.07/5.08/5.13 秒，优化后 2.86/2.87/2.88 秒。
四线程的整体进程峰值内存略有增加；搜索带紧凑存储不代表所有并发
场景的总 RSS 都下降。上述改善只代表本机和该样本，不外推其他平台。

验证包括 SIMD 与标量实现的 48 组逐单元回溯指针对照（长度 1–4096、
包含边界和整数回绕）、600 个原版粗比对样例以及精比对/BAQ 回放。
优化前后的 40,000 条完整 SAM 记录及标签一致，1/4 线程输出也一致；
检查器忽略运行命令等元数据头部。与 Stampy 的剩余差异维持上表数值，
本次加速没有提升或降低既有还原程度。哈希与种子查找代码未修改。

## 许可证

GPLv3。完整条款见 [LICENSE](LICENSE)。
