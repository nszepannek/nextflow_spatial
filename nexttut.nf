process SpaceRanger {
    tag "Run ${sample_id}" // Tagname zur Übersichtlichkeit

    input:
    val sample_id
    path transcriptome
    path probe_set
    path fastq_dir
    path image_file

    output:
    path "${sample_id}/outs", emit: spaceranger_results // Name zum weiterverwenden im workflow

    script:
    """
    ~/Masterarbeit/Spaceranger/spaceranger-3.1.3/spaceranger count \\
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
    path "umap_coordinates.csv"
    path "plots_qc"

    script:
    """
    Rscript ${r_script} ${seurat_dir}
    """
}

process Plotting_Clusters {
    tag "Plotting from ${sample_id}"

    input:
    val sample_id
    path seurat_dir2
    path r_script

    output:
    path "plots"

    script:
    """
    Rscript ${r_script} ${seurat_dir2}
    """
}

process Annotate_Data {
    tag "Annotation for ${sample_id}"

    input:
    val sample_id
    path seurat_file          // seurat_obj_with_umap.rds
    path reference_file       // .rds oder .h5ad
    path r_script

    output:
    path "csv", emit: csv_output
    path "plots_annotation", emit: annotation_plots

    when:
    reference_file.name.endsWith('.rds') || reference_file.name.endsWith('.h5ad')

    script:
    """
    Rscript ${r_script} ${seurat_file} ${reference_file}
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
    
    annot_input_ch = params.annotation_ref ? Channel.value(file(params.annotation_ref)) : Channel.empty()
    
    annot_script_ch = Channel.value(file("scripts/Annotation_Plots.R"))


    spaceranger_results = SpaceRanger(
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
    
    / Optional: Annotation only when input file exists
    Annotate_Data(
        sample_id_ch,
        clustering_results.seurat_umap,
        annot_ref_ch,
        annot_script_ch
    )

}

