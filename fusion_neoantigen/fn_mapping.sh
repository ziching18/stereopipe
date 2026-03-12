cd /stereoseq/

ulimit -n 10240
ulimit -v 33170449147
NUMBA_CACHE_DIR=/stereoseq/tmp

dataDir=/stereoseq/SD507_SD693
tmpDir=/stereoseq/tmp
refDir=/stereoseq/reference

export SINGULARITY_BIND=$dataDir,$refDir

bash /stereoseq/code/stereopipe/fusion_neoantigen/fn_saw_mapping.sh \
    -sif ~/SAW_6.0.sif \
    -genomeSize 5 \
    -splitCount 1 \
    -maskFile $dataDir/mask/SS200000873BR_D5.barcodeToPos.h5 \
    -fq1 $dataDir/reads/E100052965_r1.fq.gz \
    -fq2 $dataDir/reads/E100052965_r2.fq.gz \
    -speciesName homo_sapiens \
    -tissueType breast_tumour \
    -refIndex $refDir/STAR_SJ100 \
    -annotationFile $refDir/GRCh38.109.gtf \
    -threads 16 \
    --outTmpDir $tmpDir \
    -outDir $dataDir/result 