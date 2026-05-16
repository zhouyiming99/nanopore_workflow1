nanopolish index -d /data1/zhouyiming/qian_3/WHXWZB-2023080124A/raw_data/Nanopore/T-2/20230818-NPL231153-P4-PAQ94845/PAQ94845/fast5_pass /data1/zhouyiming/qian_3/WHXWZB-2023080124A/raw_data/Nanopore/T-2/20230818-NPL231153-P4-PAQ94845/PAQ94845/T-2.fastq
nanopolish eventalign --reads /data1/zhouyiming/qian_3/WHXWZB-2023080124A/raw_data/Nanopore/T-1/20230818-NPL231152-P4-PAQ94711/PAQ94711/T-1.fastq \
--bam /data1/zhouyiming/qian_3/drsseqres/transcript/fq/T-1_fq/T-1.bam \
--genome /data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/mRNA/osa_IRGSP_1.txdb.fa \
--scale-events \
--summary T-1_summary.txt \
--signal-index \
--threads 50 > T-1_eventalign.txt

m6anet dataprep --eventalign /data1/zhouyiming/DT4/DT4-2-1/DT4-2_eventalign.txt \
--out_dir dataprep \
--n_processes 8 \
--readcount_max 2000000

m6anet inference --input_dir dataprep --out_dir run  --pretrained_model arabidopsis_RNA002 --n_processes 8 --num_iterations 1000
