#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*  
    Usage:
	nextflow run main.nf 
    Available parameter for demuxlet v2:
	--bam <bam_file> 
	--vcf <vcf_file>
	--field <GT or GP or PL>
	--demux_tool <demuxlet or souporcell>
        --outdir <path>

container "${ params.demux_tool == 'demuxlet' ? 'popscle' : 'sourporcell' }"

*/


process demuxlet{

    publishDir "$params.outdir/demux", mode: 'copy'
    container 'conda'

    input:
	path input_bam
	path input_vcf
	val input_field

    output:
	file 'demux_out*'

    script:
    
    """
    popscle demuxlet --sam ${input_bam} --vcf ${input_vcf} --field ${input_field} --out demux_out
    """  
}



process ana_demuxlet{

    publishDir "$params.outdir/demux", mode: 'copy'
    echo true

    input:
	file demux_out
    output:
	file "violinplot.png"
	stdout

    script:
    """
    Rscript $projectDir/analysis_demuxlet.R $demux_out
    """



}

workflow demultiplex{

    if (params.demux_tool == "demuxlet"){
	    input_bam = Channel.fromPath(params.bam)
	    input_vcf = Channel.fromPath(params.vcf)
	    input_field = channel.value(params.field)

	    demuxlet(input_bam,input_vcf,input_field)
	    ana_demuxlet(demuxlet.out)
	
	    
	    
    	    
   }
  
   else if (params.demux_tool == "sourporcell"){





}
   
    	
}

workflow{
	demultiplex()
}




