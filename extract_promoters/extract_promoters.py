#!/usr/bin/env python3
# 需要安装biopython
# pip install biopython
import sys
from Bio import SeqIO

def main():
    argc = len(sys.argv)
    if argc < 5 or argc > 6:
        print(f"Usage: python {sys.argv[0]} <genome_fasta> <gff_file> <gene_id_file> <output_fa> [promoter_length]")
        sys.exit(1)

    genome_fa = sys.argv[1]    # 基因组FASTA文件
    gff_file = sys.argv[2]     # GFF3注释文件
    gene_file = sys.argv[3]    # 目标基因ID文件
    out_fa = sys.argv[4]       # 输出启动子序列文件

    # 可选启动子长度参数，默认2000bp
    prom_len = int(sys.argv[5]) if argc == 6 else 2000

    # 读取目标基因ID
    targets = set()
    with open(gene_file) as f:
        for line in f:
            tid = line.strip()
            if tid:
                targets.add(tid)

    # 载入基因组序列
    genome = SeqIO.to_dict(SeqIO.parse(genome_fa, "fasta"))

    # 解析GFF并提取启动子
    with open(out_fa, "w") as out, open(gff_file) as gff:
        for l in gff:
            if l.startswith("#"):
                continue
            c = l.rstrip().split("\t")
            if len(c) < 9 or c[2] != "mRNA":
                continue
            attrs = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
            gid = attrs.get("ID", "")
            if gid not in targets:
                continue

            chrom, start, end, strand = c[0], int(c[3]), int(c[4]), c[6]

            if strand == "+":
                s, e = max(0, start - prom_len - 1), start - 1
            else:
                s, e = end, end + prom_len

            seq = genome[chrom].seq[s:e]
            if strand == "-":
                seq = seq.reverse_complement()

            out.write(f">{gid}|{chrom}:{s+1}-{e}({strand})\n{seq}\n")

    print(f"Done, results saved to {out_fa}. Promoter length used: {prom_len} bp")

if __name__ == "__main__":
    main()
