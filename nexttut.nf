process SpaceRanger {
    spaceranger count --id="Visium_FFPE_Mouse_Brain"           
        --transcriptome=refdata-gex-mm10-2020-A           
        --probe-set=Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv           
        --fastqs=datasets/Visium_FFPE_Mouse_BraiAn_fastqs           
        --image=datasets/Visium_FFPE_Mouse_Brain_image.jpg           
        --slide=V11J26-127           
        --area=B1           
        --reorient-images=true           
        --localcores=16           
        --localmem=128
}

workflow {
    SpaceRanger
}
