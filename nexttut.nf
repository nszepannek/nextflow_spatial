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
    path spaceranger_results
    path r_script

    output:
    path "seurat_obj_with_umap.rds"
    path "umap_coordinates.csv"
    path "results/plots", optional: true

    script:
    """
    mkdir -p results/plots
    mkdir -p results/csv
    Rscript ${params.rscript_cluster} ${spaceranger_results}
    """
}

process PlotSeurat {

    tag "Plotting ${sample_id}"

    input:
    val sample_id
    path seurat_rds
    path r_script

    output:
    path "results/plots"
    path "results/csv", optional: true

    script:
    """
    mkdir -p results/plots
    mkdir -p results/csv
    Rscript ${r_script} ${seurat_rds}
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


    spaceranger_results = SpaceRanger(
        sample_id_ch,
        transcriptome_ch,
        probe_set_ch,
        fastq_dir_ch,
        image_file_ch
    )
    
    Clustering_analysis(
        sample_id_ch,
        spaceranger_results,
        seurat_script_ch
    )
    
    PlotSeurat(
    sample_id_ch,
    seurat_outputs,
    seurat_script2_ch
    )


}

