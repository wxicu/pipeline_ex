#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*  
    Usage:
	nextflow run main.nf 

    Global parameter:
	--outdir <path>
	--demux_tool <demuxlet or souporcell or vireo>
    Available parameter for demuxlet v2:
	--bam <bam_file>
	--vcf <vcf_file>
	--field <GT or GP or PL>
	
    Available paramter for souporcell:
	--bam <bam_file>
	--barcode <barcode tsv>
        --ref <reference fasta>
        --t <thread>
	--k <cluster>

   Available paramter for vireo:
	--celldata <celldata vcf>
	--ndonor <number of donor>
        --donorfile <donor vcf>
        --forceLearnGT <True or False>
	
  

*/


process demuxlet{
    publishDir "$params.outdir/demux", mode: 'copy'
    container 'popscle'

    input:
	path input_bam
	path input_vcf
	val input_field

    output:
	file 'demux*'

    script:
    
    """
    popscle demuxlet --sam ${input_bam} --vcf ${input_vcf} --field ${input_field} --out demux
    """  
}

process souporcell{

    publishDir "$params.outdir/souporcell", mode: 'copy'
    echo true

    input:
	path input_bam
	path input_barcode
	path input_ref
	val input_thread
	val input_cluster

    output:
	stdout

    script:
    
    """
    singularity exec souporcell_latest.sif souporcell_pipeline.py -i ${input_bam} -b ${input_barcode} -f ${input_ref} -t ${input_thread} -o souporcell -k ${input_cluster}
 
    """  
}


process vireo{
    publishDir "$params.outdir/vireo", mode: 'copy'
    echo true 
    /*
    container "vireo"
    */

    input:
	path input_celldata
	val input_ndonor
	path input_donorfile
	val input_forceLearnGT

    output:
	path "vireo*"
	
    

    script:
    def donor = input_donorfile.name != 'NO_FILE' ? "-d $input_donorfile" : ''
    def out = input_donorfile.name != 'NO_FILE' ? "noGT" : 'withGT'
    def learnGT = input_forceLearnGT != 'False' ? "--forceLearnGT" : ''

    """
    vireo -c ${input_celldata} -N ${input_ndonor} -o vireo_${out} ${donor} ${learnGT}
    """




}

process ana_demuxlet{

    publishDir "$params.outdir/demux", mode: 'copy'
    echo true

    input:
	file demux_out
    output:
	file "violinplot_demux.png"
	stdout

    script:
    """
    analysis_demuxlet.R -t demuxlet $demux_out
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
  
   else if (params.demux_tool == "souporcell"){
           input_bam = Channel.fromPath(params.bam)
	   input_barcode = Channel.fromPath(params.barcode)
	   input_ref = Channel.fromPath(params.ref)
	   input_thread = channel.value(params.t)
	   input_cluster = channel.value(params.k)
 	   
	   souporcell(input_bam, input_barcode, input_ref, input_thread, input_cluster)	 

}
   else if (params.demux_tool == "vireo"){
	   input_celldata = Channel.fromPath(params.celldata)
	   input_ndonor = channel.value(params.ndonor)
	   input_donorfile = Channel.fromPath(params.donorfile)
	   input_forceLearnGT = channel.value(params.forceLearnGT)
	   vireo(input_celldata,input_ndonor,input_donorfile,input_forceLearnGT)


}
    	
}

workflow{
	demultiplex()
}




