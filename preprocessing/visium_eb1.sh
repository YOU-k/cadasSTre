cd /home/yueyou/data/fastq/10x
/home/yueyou/software/spaceranger-2.1.0/spaceranger count --id="embryo1_output" \
                   --description="Adult Mouse Embryo1" \
                   --transcriptome=./refdata-gex-mm10-2020-A \
                   --fastqs=./embryo1 \
                   --image=./images/WK-XHDX-C.tif \
                   --slide=V19L01-041 \
                   --area=C1 \
                   --localcores=16 \
                   --localmem=128
