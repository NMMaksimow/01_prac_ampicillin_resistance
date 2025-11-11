#!/usr/bin/env bash
set -euo pipefail

# === PATHS ===
ALIGN_DIR="results/bwa_alignments"
STATS_DIR="results/alignment_stats"
LOG_DIR="logs"
OUTFILE="${STATS_DIR}/alignment_summary.tsv"

mkdir -p "$STATS_DIR" "$LOG_DIR"

echo -e "sample\ttotal_reads\tmapped_reads\tproperly_paired\tpercent_mapped\tpercent_proper\tmean_depth\tavg_MAPQ" > "$OUTFILE"

# === MAIN LOOP ===
for bam in ${ALIGN_DIR}/*_sorted.bam; do
    base=$(basename "$bam" _sorted.bam)
    flagfile="${STATS_DIR}/${base}_flagstat.txt"
    idxfile="${STATS_DIR}/${base}_idxstats.txt"

    echo "Processing ${base}..."

    # total reads
    total=$(grep "in total" "$flagfile" | awk '{print $1}')
    # mapped reads
    mapped=$(grep " mapped (" "$flagfile" | head -n1 | awk '{print $1}')
    # properly paired reads
    proper=$(grep "properly paired" "$flagfile" | awk '{print $1}')
    # % mapped
    pmapped=$(awk -v m=$mapped -v t=$total 'BEGIN{printf "%.2f", (m/t)*100}')
    # % proper
    pproper=$(awk -v p=$proper -v t=$total 'BEGIN{printf "%.2f", (p/t)*100}')

    # mean depth (все позиции, включая нулевые)
    depth=$(samtools depth -a "$bam" | awk '{sum+=$3; n++} END{if(n>0) print sum/n; else print 0}')

    # среднее значение MAPQ
    mapq=$(samtools view "$bam" | awk '{m[$5]++} END{for (k in m){sum+=k*m[k]; n+=m[k]} if(n>0) printf "%.2f", sum/n; else print 0}')

    echo -e "${base}\t${total}\t${mapped}\t${proper}\t${pmapped}\t${pproper}\t${depth}\t${mapq}" >> "$OUTFILE"
done

echo "Summary written to: $OUTFILE"
