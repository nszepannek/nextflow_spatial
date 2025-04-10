// Default parameter input
params.str = "Hello world!"

// splitString process
process splitString {
    publishDir "results/lower"
    
    input:
    val x

    publishDir 'results', mode: 'copy'
    
    output:
    path 'chunk_*'

    script:
    """
    printf '${x}' | split -b 6 - chunk_
    """
}

process convertToUpper {
    publishDir "results/upper"
    tag "$y"

    input:
    path y

    publishDir 'results', mode: 'copy' 

    output:
    path 'upper_*'

    script:
    """
    rev $y > upper_${y}
    """
}

// Workflow block
workflow {
    ch_str = Channel.of(params.str)     // Create a channel using parameter input
    ch_chunks = splitString(ch_str)     // Split string into chunks and create a named channel
    convertToUpper(ch_chunks.flatten()) // Convert lowercase letters to uppercase letters
}
