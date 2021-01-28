#!/bin/sh
#BSUB -q priority
#BSUB -J bcbio_flanagan
#BSUB -oo bcbio_run2.out
#BSUB -n 1
#BSUB -R "rusage[mem=8024]"
#BSUB -W 100:00
#BSUB -N rkhetani@hsph.harvard.edu

cd /n/data1/cores/bcbio/flanagan-chipseq2rna/rnaseq/flanagan_rnaseq/work/

export PATH=/opt/bcbio/local/bin:$PATH

bcbio_nextgen.py ../config/flanagan_rnaseq.yaml \
-r mincores=2 -r minconcores=2 -n 64 -t ipython -s lsf -q mcore '-rW=72:00' --retries 3 --timeout 380
