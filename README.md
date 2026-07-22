# probebwa

[![CI](https://github.com/GUIBA-EX/probebwa/actions/workflows/ci.yml/badge.svg)](https://github.com/GUIBA-EX/probebwa/actions/workflows/ci.yml)

UCE(超保守元件)探针/bait 序列比对到参考基因组的短读比对器:密集 k-mer 种子
扩展、banded affine-gap 比对(质量感知错配罚分 + 同聚物感知 gap-open)、贝叶斯
后验 MAPQ 打分;支持单端和配对末端(候选 shortlist + 联合后验,能把落进重复
区域的读段纠正到正确拷贝)。

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

## 许可证

GPLv3。完整条款见 [LICENSE](LICENSE)。
