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
    ${params.spaceranger_path} count \\
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
    path "${params.outdir}/plots_qc/${sample_id}", mode: 'copy'

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
    path "${params.outdir}/plots/${sample_id}", mode: 'copy'

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
    path reference_file       // .rds oder .h5ad
    path r_script

    output:
    path "${params.outdir}/csv/${sample_id}", mode: 'copy', emit: csv_output
    path "${params.outdir}/plots_annotation/${sample_id}", mode: 'copy', emit: annotation_plots

    when:
    reference_file.name.endsWith('.rds') || reference_file.name.endsWith('.h5ad')

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
    
    
    annot_ref_ch = params.annot_ref ? Channel.value(file(params.annot_ref)) : Channel.empty()
    
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
    
    
   
    Annotate_Data(
        sample_id_ch,
        clustering_results.seurat_umap,
        annot_ref_ch,
        annot_script_ch
    )

}

