#!/usr/bin/env bash

mkdir -p input/datasets

cd input

# transcript_ref
curl -O https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz
tar -xzvf refdata-gex-mm10-2020-A.tar.gz && rm refdata-gex-mm10-2020-A.tar.gz
mv refdata-gex-mm10-2020-A transcript_ref

# probeset
curl https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv -o probe_set.csv

# datasets/fastq_data
curl https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Mouse_Brain/Visium_FFPE_Mouse_Brain_fastqs.tar -o datasets/Visium_FFPE_Mouse_Brain_fastqs.tar
tar -xvf datasets/Visium_FFPE_Mouse_Brain_fastqs.tar -C datasets/ && rm datasets/Visium_FFPE_Mouse_Brain_fastqs.tar
mv datasets/Visium_FFPE_Mouse_Brain_fastqs datasets/fastq_data

# datasets/input_image.jpg
curl https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Mouse_Brain/Visium_FFPE_Mouse_Brain_image.jpg -o datasets/input_image.jpg

# annot_ref.h5ad
curl https://datasets.cellxgene.cziscience.com/f87d516e-83fa-4ca4-a37a-ed1ab7ff2199.h5ad -o annot_ref.h5ad