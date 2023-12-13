#!/bin/bash

# 初始化参数变量
genome=""
gff3=""
cds=""
pep=""
id=""
up=2000
down=2000
compress=false  # Flag to indicate whether to compress

# 显示帮助信息的函数
show_help() {
	echo "Usage: $0 [--genome genome_file] [--gff3 gff3_file] [--cds cds_file] [--pep pep_file] [--id gene/mRNA_id] [--up up] [--down down] [--gz]"
	echo "Options:"
	echo "  --genome  Specify the genome fasta file"
	echo "  --gff3    Specify the gff3 file"
	echo "  --cds     Specify the cds fasta file"
	echo "  --pep     Specify the pep fasta file"
	echo "  --id      Specify the gene/mRNA id"
	echo "  --up      Specify how many bp upstream for gene/mRNA"
	echo "  --down    Specify how many bp downstream for gene/mRNA"
	echo "  --gz      Compress all the result if this option is present"
	echo "  --help    Display this help message"
}

# Check if no arguments are provided
if [ $# -eq 0 ]; then
	show_help
	exit 0
fi

# 循环处理参数
while [[ $# -gt 0 ]]; do
	case "$1" in
		--genome)
			shift
			genome="$1"
			;;
		--gff3)
			shift
			gff3="$1"
			;;
		--cds)
			shift
			cds="$1"
			;;
		--pep)
			shift
			pep="$1"
			;;
		--id)
			shift
			id="$1"
			;;
		--up)
			shift
			if [ -n "$1" ]; then
				up="$1"
			else
				echo "Error: --up requires a non-empty value."
				exit 1
			fi
			;;
		--down)
			shift
			if [ -n "$1" ]; then
				down="$1"
			else
				echo "Error: --down requires a non-empty value."
				exit 1
			fi
			;;
		--gz)
			compress=true
			;;
		--help)
			show_help
			exit 0
			;;
		*)
			echo "Unknown option: $1"
			show_help
			exit 1
			;;
	esac
	shift
done

grep -E "ID=$id;|ID=$id$|Parent=$id;|Parent=$id$" $gff3 > $id.gff
line=$(grep -E "ID=$id;|ID=$id$" $gff3)

# 使用awk来提取列并进行计算
chr=$(echo "$line" | awk '{print $1}')
start=$(echo "$line" | awk '{print $4 - $up}')
end=$(echo "$line" | awk '{print $5 + $down}')

# 提取gff3
grep -E "ID=$id;|ID=$id$|Parent=$id;|Parent=$id$" $gff3 > $id.gff
# 提取基因组序列
if [ -e "$genome.fai" ]; then
	samtools faidx $genome $chr:$start-$end > $id.genome_plus_up${up}_down${down}.fa
else
	samtools faidx $genome
	samtools faidx $genome $chr:$start-$end > $id.genome_plus_up${up}_down${down}.fa
fi
# 提取cds
if [ -e "$cds.fai" ]; then
	samtools faidx $cds $id > $id.cds.fa
else
	samtools faidx $cds
	samtools faidx $cds $id > $id.cds.fa
fi
# 提取pep
if [ -e "$cds.fai" ]; then
	samtools faidx $pep $id > $id.pep.fa
else
	samtools faidx $pep
	samtools faidx $pep $id > $id.pep.fa
fi

# Check if --gz is present, then compress all result
if [ "$compress" = true ]; then
	echo "Compressing all the result..."
	tar zcvf $id.tar.gz $id.genome_plus_up${up}_down${down}.fa $id.gff $id.cds.fa $id.pep.fa
	rm $id.genome_plus_up${up}_down${down}.fa $id.gff $id.cds.fa $id.pep.fa
else
	echo "Not compressing the result"
fi

