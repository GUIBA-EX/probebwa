# probebwa

一个用 Rust 写的短读比对器,专门用于 UCE(超保守元件)探针/bait 设计流程里的
比对环节——把候选探针/bait 序列比对到参考基因组,评估唯一性、检查是否落在
重复区域、确认候选位点。核心是密集 k-mer 种子扩展、banded affine-gap 比对、
贝叶斯后验 MAPQ 打分。

已经在真实数据上验证过(见[真实数据验证](#真实数据验证)),但还没在染色体级
基因组(人类等)上测过;见[已知局限](#已知局限)。

## 算法

1. **种子**——读段里每个重叠的 15-mer 都精确查找;单点错配邻居按读长分档,
   只在一部分位置补充探测,再经过一层轻量的碱基组成相似度过滤,才作为候选保留。
   含 N(或其他非 ACGT 模糊码)的窗口直接跳过,不参与查找——参考基因组建索引时
   同样跳过含 N 的窗口,两边保持一致,否则 N 会被静默当成 'A' 处理,污染哈希表。
2. **聚类**——种子命中按链方向感知的对角线分组,用滑动窗口合并成候选位点。
3. **比对**——banded、affine-gap、质量感知 + 同聚物感知的动态规划,为每个候选
   产出 CIGAR 和得分。
4. **打分**——候选之间的贝叶斯后验(按每个碱基质量、沿 CIGAR 走一遍算出)结合
   漏找位点概率和随机匹配显著性检验,得到 MAPQ。
5. **配对末端**——每条 read 先各自建候选位点列表和后验概率,再按覆盖 99.9% 单端
   后验质量的候选(3~20 个)建 shortlist,把配偶直接比对到插入片段大小推出的窗口上,
   对所有候选对算联合似然(两端似然 + 插入片段大小的高斯项,异常间距/方向会被
   一个可配置的结构变异先验兜底),挑最好的一对。这样即使一条 read 自己的最佳
   候选选错了(比如落在重复序列的另一个拷贝上),只要配偶能唯一锚定,还是能被
   纠正到正确的拷贝上。proper pair 要求 FR 方向且间距落在插入片段大小估计的
   5 个标准差以内(该估计在运行中在线更新)。插入片段大小统一按 SAM TLEN 语义
   计算(上游 read 最左端到下游 read 最右端的完整跨度,经 CIGAR 换算,而不是两端
   起点的简单差值),候选打分、mate rescue 窗口定位、在线学习用的都是同一个口径。

这几块模型(种子扫描的长度自适应采样、漏找位点概率、随机匹配显著性检验)都是
针对探针/bait 唯一性评估这个具体场景独立设计实现的,不基于任何专有源码。
哪里做了简化,见[已知局限](#已知局限)。

## 编译

```bash
cargo build --release
```

默认按编译机器实际支持的 CPU 特性编译(`.cargo/config.toml` 里的
`-C target-cpu=native`,由 rustc/LLVM 在编译时检测,不需要额外的检测逻辑)——
有 AVX-512 之类的指令集时,编译器能更激进地自动向量化基因组解包、碱基计数这类
热点循环。代价是产出的二进制跟编译机器的 CPU 绑定,拿到特性更少的机器上跑会
直接因为非法指令崩溃,而不是默默降级变慢。如果需要编译一次、拷贝到别的机器上
跑,用更保守的基线覆盖掉默认值:

```bash
RUSTFLAGS="-C target-cpu=x86-64-v2" cargo build --release
```

(或目标机器实际能保证的任何基线特性集)。

## 用法

```bash
# 1. 建基因组索引(.stidx)
./target/release/probebwa build-genome \
    --species human --assembly hg38 \
    -G hg38 ref.fa.gz

# 2. 建 15-mer 哈希表(.sthash)
./target/release/probebwa build-hash --genome hg38 -H hg38

# 3a. 比对探针/bait 序列(单端)
./target/release/probebwa map --genome hg38 --hash hg38 \
    --substitution-rate 0.001 -M probes.fq.gz > output.sam

# 3b. 比对配对末端 reads(两个文件即视为配对)
./target/release/probebwa map --genome hg38 --hash hg38 \
    -M reads_1.fq.gz reads_2.fq.gz > output.sam

# 其他选项
./target/release/probebwa map --genome hg38 --hash hg38 --phred64 -M old.fq > out.sam
./target/release/probebwa map --genome hg38 --hash hg38 --inputformat=fasta -M probes.fa > out.sam
./target/release/probebwa map --genome hg38 --hash hg38 \
    --readgroup "ID:rg1,SM:sample1,PL:illumina" -M probes.fq.gz > out.sam

# 多线程(rayon 线程池;0 = 用满所有核心,默认 1 即严格串行)
./target/release/probebwa map --genome hg38 --hash hg38 \
    --threads 8 -M reads_1.fq.gz reads_2.fq.gz > out.sam

# gap 罚分(Phred 标度,默认 --gapopen 40 --gapextend 3,内部按本项目自己校准过
# 的基线换算,细节见 align::smith_waterman 的模块文档)
./target/release/probebwa map --genome hg38 --hash hg38 \
    --gapopen 30 --gapextend 2 -M probes.fq.gz > out.sam

# 调整结构变异先验(配对末端异常间距/方向的兜底概率,Phred 标度,默认 55)
./target/release/probebwa map --genome hg38 --hash hg38 \
    --svprior 45 -M reads_1.fq.gz reads_2.fq.gz > out.sam

# 接受旧版(pre-CASAVA-1.8)`/1`/`/2` 后缀的 read ID(现代 CASAVA 1.8+ 头部
# 不受此开关影响,始终按空格分字段正确解析)
./target/release/probebwa map --genome hg38 --hash hg38 \
    --casava8 -M old_reads_1.fq old_reads_2.fq > out.sam

# 比对前从每条读段 3' 端剪掉接头序列(容忍测序错误:每 10bp 重叠允许 1 个错配)
./target/release/probebwa map --genome hg38 --hash hg38 \
    --adapter-strip AGATCGGAAGAGC -M probes.fq.gz > out.sam

# BAM 输出(必须指定 --output;二进制不该直接输出到终端)
./target/release/probebwa map --genome hg38 --hash hg38 \
    --outputformat bam --output out.bam -M probes.fq.gz

# BAM + 坐标排序 + 生成 .bai 索引(要求 --outputformat=bam)
./target/release/probebwa map --genome hg38 --hash hg38 \
    --outputformat bam --output out.bam --index -M probes.fq.gz

# 只保留 read ID 以指定前缀开头的记录
./target/release/probebwa map --genome hg38 --hash hg38 \
    --labelfilter sample1_ -M probes.fq.gz > out.sam
```

FASTA 输入(`--inputformat=fasta`)是 UCE 探针/bait 集合的常见形式,这类数据通常
以未比对的 FASTA 分发;探针会拿到一条统一的高置信度质量值,因为本来就没有真实的
测序质量数据。

## 性能

- **多线程**——单端走 `rayon` 的 `par_iter`(`--threads` 控制池大小,0 = 用满所有
  核心);配对末端因为插入片段大小分布是在线学习、天然串行,所以分批处理:
  开头一小段串行跑起来"预热"模型(避免第一批完全拿着默认先验打分),之后按批
  冻结模型快照、批内并行,批与批之间再把观测结果串行折算回模型(`src/lib.rs`
  的 `PAIR_WARMUP_SIZE`/`PAIR_BATCH_SIZE` 有完整推导);每一对配偶自己的种子+
  聚类+比对这一步,两条 read 也用 `rayon::join` 并行。
- **`mimalloc` 全局分配器**——DP 比对和候选列表这些热路径每次都会有不少小块、
  短生命周期的 `Vec` 分配,`mimalloc` 在这种模式下比系统分配器稳定更快
  (`src/main.rs`)。
- **编译期 CPU 特性自适应**——见上面[编译](#编译)一节的 `target-cpu=native`。
- **`profile.release`**——`lto = true`、`codegen-units = 1`、`panic = "abort"`
  (`Cargo.toml`),牺牲编译时间换运行时性能和更小的二进制。

这些改动全部是"语义不变、只改编译器/运行时决策"的性能优化,每一步都拿
`cargo test --release` + 真实 E. coli 数据的逐字节 diff 校验过没有引入行为差异
(唯一一次尝试性改动——DP 缓冲区跨调用复用——测出来反而更慢,已经回退,细节见
`align::smith_waterman` 的模块文档)。在本项目自己的 2 万对 E. coli 真实数据
基准上,这一整套优化把单机墙钟时间从约 58 秒降到约 1.3 秒。

## 架构

| 模块 | 作用 |
|---|---|
| `index/` | 基因组索引(`.stidx`,2-bit 打包)+ 15-mer 哈希表(`.sthash`) |
| `mapper/` | 种子扫描、候选聚类、单端/配对末端比对 |
| `align/` | Banded affine-gap 比对、CIGAR 工具(解析/格式化/参考跨度计算) |
| `mapq/` | 贝叶斯后验、漏找位点概率、MAPQ、插入片段大小在线学习 |
| `io/` | FASTA/FASTQ 解析(含 `--casava8`/`--adapter-strip` 预处理)、SAM/BAM 输出 |

## 设计要点

延续 [minimap2](https://github.com/lh3/minimap2) 这类生产级比对器的一般做法:

- **2-bit 打包基因组**(`types::Contig`)——4 个碱基压成 1 字节,N 游程单独存储,
  解码时再还原。
- **零分配种子探测**(`mapper::seeds`)——错配邻居生成靠原地改写一块栈缓冲区,
  而不是每个变体都分配一次 `Vec<u8>`。
- **复用 DP 缓冲区**(`align::smith_waterman`)——行缓冲区互相交换复用,traceback
  矩阵是一整块扁平缓冲区,不是每行一个 `Vec`。
- **重复感知的种子探测**——一个位置的精确匹配命中数已经看起来很重复时,就跳过
  错配变体探测。

## 真实数据验证

用真实 E. coli(K-12 MG1655)参考基因组 + 真实 ENA 测序数据(DRR002055,2 万条
配对末端 read)做过验证,并和另一个成熟短读比对工具在同一份数据上的输出做过
逐条交叉比较,而不是只靠合成测试。这套流程在开发过程中实际抓出过几个真实 bug
(比如 gap-open 罚分基线偏低导致在重复样区域用虚假 indel 换掉真实错配、TLEN
计算漏算下游 read 自身的参考跨度、配对候选打分和插入片段大小模型的口径不一致、
参考基因组里的 N 被静默当 'A' 索引),都是靠这套真实数据交叉验证才发现的,不是
单纯的合成数据能覆盖到的。目前在这份基准数据上,和对照工具的比对位点一致率
维持在高位,两边都给出高置信度却选了不同位点的情况是少数。

这套验证不是自动化测试的一部分(需要真实参考基因组和外部对照工具,不适合放进
CI),但每一次涉及打分/评分逻辑的改动都建议照这个流程重新跑一遍。

## 测试

```bash
cargo test --release
```

- `src/index/hashtable.rs` 里的单元测试——N-gap 相关的哈希表建库/查询回归测试
  (`n_windows_are_not_indexed`、`is_acgt_only_rejects_n_and_other_ambiguity_codes`)。
- `tests/integration.rs`——核心数据结构的单元测试(k-mer 哈希、反向互补、质量解析)。
- `tests/mapping.rs`——走 CLI 同一套公开 API 的端到端测试:链方向定位、错配、
  indel CIGAR、未比对读段、proper/improper 配对、mate rescue、N 游程,以及一个
  300 条读段/20kb 的合成准确率检查(±2bp 内 100% 正确;加 `-- --nocapture` 能打印
  出具体数字)。这部分都是合成数据,真实数据的验证见上一节。
- `tests/cli.rs`——针对只存在于 `main.rs`(而非库)里的开关(`--outputformat=bam`、
  `--labelfilter`)运行编译出的二进制做验证。
- `examples/tune_homopolymer.rs`——一个手动跑的调参脚手架,不算在测试套件里
  (见下方已知局限)。
- `examples/profile_paired.rs`——一个手动跑的计时脚手架,拿真实 FASTQ 数据分别
  给配对末端映射的"共享种子+比对"部分和"shortlist/交叉比对"部分计时,定位
  性能改动该往哪投入;同样不算在测试套件里,用法见文件头注释。

## 已知局限

- **还没在染色体级基因组(人类等)上测过。** 真实数据验证目前只覆盖 E. coli
  这个量级的基因组;人类基因组规模下哈希表内存占用、候选数量、`--threads`
  的扩展性都还没有实测数据。
- **MAPQ 打分模型没有折算逐位置的重复掩码,长读段用的是简化公式,不是完整的
  子区间最小化处理。** 一个能看到的具体影响是:一条完美、唯一比对的 ~40bp 读段
  报出的 MAPQ 是 ≈80 而不是固定的 99——因为短读段每个相位只有 2 个独立候选窗口,
  这是符合模型本身逻辑的行为,不是退步。
- **`align::smith_waterman` 是标量实现,没上 SIMD——而且是故意不上。** 现成的
  vetted 库方案(`block-aligner`)只支持固定的匹配/错配和 gap 罚分;而这个比对器
  质量感知的错配罚分、同聚物感知的 gap-open 正是探针唯一性评估这个场景里值得
  保留的行为,所以没有为了 SIMD 吞吐量牺牲它们。真要在不退步的前提下拿到 SIMD,
  得自己针对 block-aligner 的 unsafe 内部结构写一个核酸专用的打分 profile——
  工作量不小,而且某种意义上还是"自己手撸 SIMD"的那类风险,只是换了一层。
- **DP 里有两个调参常数是经验搜索出来的。** `align::smith_waterman` 里同聚物
  gap-open 折扣的斜率、错配罚分质量限幅的下限,是拿一套合成的"同聚物 indel +
  噪声"测试(`examples/tune_homopolymer.rs`)扫出来的;现在用的值已经在扫描
  结果的最优平台上。
- **配对末端候选 shortlist + 联合后验已经实现,但有几处设计上需要自己拿主意
  的地方。** 胜出的候选对如果就是两端各自独立的最佳候选,后验质量该怎么算是
  清楚的;但两端都不是独立最佳时(被交叉比对出来的候选反超)该怎么处理,没有
  唯一"正确"答案。这个实现里(`mapper::paired::PairedEndMapper::mapq_for_winning_pair`
  有完整推导注释):
  - 两端独立最佳一致——各自报自己的单端 MAPQ;
  - 恰好一端独立最佳——两端都报那一端(锚定 read)的单端置信度,因为另一端的
    位置完全是靠锚定 read 才确定下来的;
  - 两端都不是独立最佳——这里选择用胜出候选对本身的联合后验分别给两端打分。
  另外,重复位点掩码没有折进配偶交叉比对的打分里(跟单端漏找位点概率模型那条
  已知局限是同一个遗留问题);候选对数量最坏情况下能到 20×20,单条模糊 read pair
  的计算量比只做单端独立打分时重不少。
- **`tests/mapping.rs::test_paired_disambiguates_repeat_copy_via_anchor`**
  是一个针对性回归测试:构造一条 read 唯一锚定、配偶自己比对反而会选错重复拷贝
  的场景,验证配对信息确实能把配偶纠正回正确拷贝——这正是这套机制存在的意义。

## 许可证

GPLv3。完整条款见 [LICENSE](LICENSE)。
