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
	--field <GT or GP or PL, default GT>
    --alpha <default 0.5>
    --alpha_list <grid search, default False> becasue --alpha 0 --alpha 0.5 is different from --alpha 0.5 TUT 
    --geno_error_offset <default 0.1>
	
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
	


not important
echo ${demux_files} > print.txt

*/


process demuxlet{
    publishDir "$params.outdir/demux", mode: 'copy'
    container 'popscle'

    input:
	file bam
	file vcf 
	val field
    val alpha_list
    each alpha
    each genoerror

    output:
	file 'demux*'

    script:
    def alpha_value = alpha_list.replaceAll(/,/, " --alpha ")
    if (alpha_list == 'False')
	"""
	popscle demuxlet --sam $bam --vcf $vcf --field $field --alpha $alpha --geno-error-offset $genoerror --out demux_${field}_${alpha}_${genoerror}
	"""  

    else   
	"""
	popscle demuxlet --sam $bam --vcf $vcf --field $field --alpha ${alpha_value} --geno-error-offset $genoerror --out demux_${field}_${alpha_list}_${genoerror}
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
	didnt understand why this command cannot recognize docker
    */

    input:
	path celldata
	val ndonor
	path donorfile
	val forceLearnGT

    output:
	path "vireo*"
	

    script:
    def donor = donorfile.name != 'NO_FILE' ? "-d $donorfile" : ''
    def out = donorfile.name != 'NO_FILE' ? "noGT" : 'withGT'
    def learnGT = forceLearnGT != 'False' ? "--forceLearnGT" : ''

    """
    vireo -c $celldata -N $ndonor -o vireo_${out} $donor $learnGT
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


process compare_parameter_demux{
    publishDir "$params.outdir/demux", mode: 'copy'
    echo true

    input:
	file demux_result
        
    output:
    file "result.tsv"
    file "result_detail.tsv"
    file "assignment.tsv"
	file "barplot.png"
	file "barplot by group.png"

    script:
    def demux_files = ""
    for(r : demux_result) {
    	demux_files = r + " " + demux_files
    }

    """
    compareParam.R ${demux_files}  
    """
}

workflow demultiplex{
    if (params.demux_tool == "demuxlet"){
	    input_bam = Channel.fromPath(params.bam)
	    input_vcf = Channel.fromPath(params.vcf)
	    input_field = channel.value(params.field)
        input_alpha_list = Channel.value(params.alpha_list)

        if (params.alpha =~ /,/){
			input_alpha = Channel.from(params.alpha)
					             .map{ it.split(',') }
					             .flatten()
	    }
        else{
			input_alpha = Channel.from(params.alpha)
	    }
	      
        if (params.geno_error_offset =~ /,/ ){
			input_genoerror = Channel.from(params.geno_error_offset)
									 .map{ return it.tokenize(',')}
									 .flatten()
	    }
        else{
			input_genoerror = Channel.from(params.geno_error_offset)
	    }

	    demuxlet(input_bam, input_vcf, input_field, input_alpha_list, input_alpha, input_genoerror)
	    
		if(params.alpha_list == 'False') {
			compare_parameter_demux(demuxlet.out.collect())
	    }
	    else{
			ana_demuxlet(demuxlet.out)
	    }
     
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




