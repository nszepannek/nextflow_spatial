process Validate_Inputs {
    tag "Validate input data ${sample_id}"

    input:
    val sample_id
    path transcriptome
    path probe_set
    path fastq_dir
    path image_file

    output:
      val 'ok'

    script:
    """
        # Check that transcriptome is a file
        if [ ! -d "$transcriptome" ]; then
          echo "Error: transcriptome does not point to a valid file: $transcriptome" >&2
          exit 1
        fi

        # Check that probe_set is a file
        if [ ! -f "$probe_set" ]; then
          echo "Error: probe_set does not point to a valid file: $probe_set" >&2
          exit 1
        fi

        # Check that fastq_dir is a directory
        if [ ! -d "$fastq_dir" ]; then
          echo "Error: fastq_dir does not point to a valid directory: $fastq_dir" >&2
          exit 1
        fi

        # Check that image_file is a file
        if [ ! -f "$image_file" ]; then
          echo "Error: image_file does not point to a valid file: $image_file" >&2
          exit 1
        fi
    """
}

process SpaceRanger {
    tag "Run ${sample_id}" // Tagname zur Übersichtlichkeit

    input:
    val validated
    val sample_id
    path transcriptome
    path probe_set
    path fastq_dir
    path image_file

    output:
    path "${sample_id}/outs", emit: spaceranger_results // Name zum weiterverwenden im workflow

    script:
    """
    spaceranger count \\
        --id=${sample_id} \\
        --transcriptome=${transcriptome} \\
        --probe-set=${probe_set} \\
        --fastqs=${fastq_dir} \\
        --image=${image_file} \\
        --slide=V11J26-127 \\
        --area=B1 \\
        --reorient-images=true \\
        --localcores=${task.cpus} \\
        --localmem=128 \\
        --create-bam=false
    """
}

process Clustering_analysis {
    tag "Clustering from ${sample_id}"

    input:
    val sample_id
    path seurat_dir
    path r_script

    output:
    path "seurat_obj_with_umap.rds", emit: seurat_umap
    path "umap_coords_with_clusters.csv"
    path "${params.outdir}/plots_qc/${sample_id}"

    script:
    """
    mkdir -p ${params.outdir}/plots_qc/${sample_id}
    Rscript ${r_script} ${seurat_dir}
    cp -r plots_qc/* ${params.outdir}/plots_qc/${sample_id}/
    """
}

process Plotting_Clusters {
    tag "Plotting from ${sample_id}"

    input:
    val sample_id
    path seurat_dir2
    path r_script

    output:
    path "${params.outdir}/plots/${sample_id}"

    script:
    """
    mkdir -p ${params.outdir}/plots/${sample_id}
    Rscript ${r_script} ${seurat_dir2}
    cp -r plots/* ${params.outdir}/plots/${sample_id}/
    """
}

process Annotate_Data {
    tag "Annotation for ${sample_id}"

    input:
    val sample_id
    path seurat_file          // seurat_obj_with_umap.rds
    val reference_file       // .rds oder .h5ad
    path r_script

    output:
    path "${params.outdir}/csv/${sample_id}", emit: csv_output
    path "${params.outdir}/plots_annotation/${sample_id}", emit: annotation_plots

    script:
    """
    mkdir -p ${params.outdir}/csv/${sample_id}
    mkdir -p ${params.outdir}/plots_annotation/${sample_id}

    Rscript ${r_script} ${seurat_file} ${reference_file}

    cp -r csv/* ${params.outdir}/csv/${sample_id}/
    cp -r plots_annotation/* ${params.outdir}/plots_annotation/${sample_id}/
    """
}

workflow {
    Channel.value(params.sample_id).set { sample_id_ch }
    Channel.value(file(params.transcriptome)).set { transcriptome_ch }
    Channel.value(file(params.probe_set)).set { probe_set_ch }
    Channel.value(file(params.fastq_dir)).set { fastq_dir_ch }
    Channel.value(file(params.image_file)).set { image_file_ch }
    Channel.value(file("scripts/Clustering_UMAP.R")).set { seurat_script_ch }
    Channel.value(file("scripts/Plots_Clusters.R")).set { seurat_script2_ch }
    Channel.value(file("scripts/Annotation_Plots.R")).set { annot_script_ch }

    validate_input_result = Validate_Inputs(
        sample_id_ch,
        transcriptome_ch,
        probe_set_ch,
        fastq_dir_ch,
        image_file_ch
    )

    spaceranger_results = SpaceRanger(
        validate_input_result,
        sample_id_ch,
        transcriptome_ch,
        probe_set_ch,
        fastq_dir_ch,
        image_file_ch
    )
    
    clustering_results = Clustering_analysis(
    sample_id_ch,
    spaceranger_results,
    seurat_script_ch
    )
    
    Plotting_Clusters(
        sample_id_ch,
        clustering_results.seurat_umap,
        seurat_script2_ch
    )
    
    if (params.run_annotation) {
      annot_ref_ch = Channel.value(params.annot_ref)()
        Annotate_Data(
            sample_id_ch,
            clustering_results.seurat_umap,
            annot_ref_ch,
            annot_script_ch
        )
    }

    else {
        log.info " Annotation step skipped (params.run_annotation = false)"
    }    
}

