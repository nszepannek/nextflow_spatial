# Input data

```plaintext
input/
├── transcript_ref/
├── datasets/
│   ├── fastq_data/
│   └── input_image.jpg
├── probe_set.csv
├── annot_ref.h5ad
└── annot_ref.rds  (optional)
```

# Build container
> docker build -t nf-spatial:latest .


# Run containerized nextflow

> ./run.sh