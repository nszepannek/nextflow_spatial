process SpaceRanger {
    tag "Run ${sample_id}"

    input:
    val sample_id
    path transcriptome from file(params.transcriptome)
    path probe_set     from file(params.probe_set)
    path fastq_dir     from file(params.fastq_dir)
    path image_file    from file(params.image_file)

    output:
    path "${sample_id}" into spaceranger_results

    cpus 8
    memory '64 GB'

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

workflow {
    Channel.value(params.sample_id).set { sample_id_ch }

    SpaceRanger(sample_id_ch)
    spaceranger_results.view { "Result: $it" }
}