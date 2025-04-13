process SpaceRanger {
    tag "Run ${sample_id}"

    input:
    val sample_id
    path transcriptome
    path probe_set
    path fastq_dir
    path image_file

    output:
    path "${sample_id}/outs", emit: spaceranger_results

    cpus 8
    memory '64 GB'

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

workflow {
    Channel.value(params.sample_id).set { sample_id_ch }
    Channel.value(file(params.transcriptome)).set { transcriptome_ch }
    Channel.value(file(params.probe_set)).set { probe_set_ch }
    Channel.value(file(params.fastq_dir)).set { fastq_dir_ch }
    Channel.value(file(params.image_file)).set { image_file_ch }

    SpaceRanger (
            sample_id_ch,
            transcriptome_ch,
            probe_set_ch,
            fastq_dir_ch,
            image_file_ch
        )    spaceranger_results.view { "Result: $it" }
}