#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*  
    Usage:
	nextflow run main.nf 
    Available parameter:
	--bam <bam_file> 
	--vcf <vcf_file>
	--field <GT or GP or PL>
	

*/


process demuxlet{
    input:
	path input_bam
	path input_vcf
	val input_field

    script:
    """
    popscle demuxlet --sam ${input_bam} --vcf ${input_vcf} --field ${input_field} --out output/demux
    """
}

workflow runDemuxlet{
    input_bam = Channel.fromPath(params.bam)
    input_vcf = Channel.fromPath(params.vcf)
    input_field = channel.value(params.field)
	
    demuxlet(input_bam,input_vcf,input_field)
    	
}

workflow{

	runDemuxlet()
}




