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
	--alpha
	--alpha_list <grid search> becasue --alpha 0 --alpha 0.5 is different from --alpha 0.5 TUT 
	--geno_error_offset 
	
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
	
Some notes:
echo ${demux_files} > print.txt
singularity exec -B $PWD/souporcell/:/souporcell $projectDir/souporcell_latest.sif ls $out 
ingularity exec -B $PWD/souporcell/:/souporcell $projectDir/souporcell_latest.sif cp -R $out .
data = Channel.fromPath(data_path, type: 'dir')	 
 
*/

def split_input(input){
	if (input =~ /,/ ){
		Channel.from(input).map{ return it.tokenize(',')}.flatten()
	}
    else{
		Channel.from(input)
	}

}
process demuxlet {
    publishDir "$params.outdir/demux", mode: 'copy'
    container 'popscle'

    input:
	file sam
    each tag_group 
    each tag_UMI
    
	each plp

	file vcf 
	each field
	each geno_error_offset
	each geno_error_coeff
	each r2_info
	each min_mac
	each min_callrate
	each sm
	each sm_list
	
	val alpha_list
	each alpha
	each doublet_prior
	each sam_verbose
	each vcf_verbose

	each cap_BQ 
	each min_BQ
	each min_MQ
	each min_TD
	each excl_flag

	each group_list
	each min_total  
	each min_umi
	each min_snp


    output:
	file 'demux*'

    script:
	
    def samfile = "--sam $sam"
    def taggroup = tag_group != 'False' ? "--tag-group ${tag_group}" : ''
    def tagUMI = tag_UMI != 'False' ? "--tag-UMI ${tag_UMI}" : ''
    def plpfile = plp != 'no_plp' ? "--plp $plp" : ''
    def vcffile = "--vcf $vcf"
    def fieldinfo = "--field $field"
    def genoerror_off = "--geno-error-offset ${geno_error_offset}"
    def genoerror_cof = "--geno-error-coeff ${geno_error_coeff}"
    def r2info = "--r2-info ${r2_info}"
    def minmac = "--min-mac ${min_mac}"
    def mincallrate = "--min-callrate ${min_callrate}"	
    def smlist = sm != 'False' ? "--sm $sm" : ''
    def sm_list_file = sm_list != 'no_sm_list' ? "--sm-list ${sm_list}" : ''
    def alpha_value = alpha_list.replaceAll(/,/, " --alpha ")
    def doubletprior = "--doublet-prior ${doublet_prior}"		
    def samverbose = "--sam-verbose ${sam_verbose}"
    def vcfverbose = "--vcf-verbose ${vcf_verbose}"
    def capBQ = "--cap-BQ ${cap_BQ}"
    def minBQ = "--min-BQ ${min_BQ}"
    def minMQ = "--min-MQ ${min_MQ}"
    def minTD = "--min-TD ${min_TD}"
    def exclflag = "--excl-flag ${excl_flag}"
    def grouplist = group_list != 'False' ? "--group-list ${group_list}" : '' 
    def mintotal = "--min-total ${min_total}"
    def minumi = "--min-umi ${min_umi}"
    def minsnp = "--min-snp ${min_snp}"

    def out = "demux+${tag_group}+${tag_UMI}+${plp}+${field}+${geno_error_offset}+${geno_error_coeff}+${r2_info}+${min_mac}+${min_callrate}+${sm}+${sm_list}+${alpha}+${doublet_prior}+${sam_verbose}+${vcf_verbose}+${cap_BQ}+${min_BQ}+${min_MQ}+${min_TD}+${excl_flag}+${group_list}+${min_total}+${min_umi}+${min_snp}"
   
    if (alpha_list == 'False')
	"""
	popscle demuxlet $samfile $taggroup $tagUMI $plpfile $vcffile $fieldinfo ${genoerror_off} ${genoerror_cof} $r2info $minmac $mincallrate $smlist ${sm_list_file} --alpha $alpha $doubletprior $samverbose $vcfverbose $capBQ $minBQ $minMQ $minTD $exclflag $grouplist $mintotal $minumi $minsnp --out $out
	"""  

    else{
    out = "demux+${tag_group}+${tagUMI}+${plp}+${field}+${geno_error_offset}+${geno_error_coeff}+${r2_info}+${min_mac}+${min_callrate}+${sm}+${sm_list}+${alpha_list}+${doublet_prior}+${sam_verbose}+${vcf_verbose}+${cap_BQ}+${min_BQ}+${min_MQ}+${min_TD}+${excl_flag}+${group_list}+${min_total}+${min_umi}+${min_snp}"	
	"""
	popscle demuxlet $sam $taggroup $tagUMI $plpfile $vcffile $fieldinfo ${genoerror_off} ${genoerror_cof} $r2info $minmac $mincallrate $smlist ${sm_list_file} --alpha ${alpha_value} $doubletprior $samverbose $vcfverbose $capBQ $minBQ $minMQ $minTD $exclflag $grouplist $mintotal $minumi $minsnp --out $out
	"""  
    }
}

process souporcell{

    publishDir "$params.outdir/souporcell", mode: 'copy'
    echo true

    input:
	val bam
	val barcodes
	val fasta
	each threads
	each clusters

  	each ploidy
  	each min_alt 
  	each min_ref
	each max_loci
	each restarts
        
	each common_variants
	each known_genotypes
	each known_genotypes_sample_names 
		
	each skip_remap
	each ignore


    output:
	path 'soup*'
	
    script:	

	def bamfile = "-i $bam"
	def barcode = "-b $barcodes"
	def fastafile = "-f $fasta"
	def thread = "-t $threads"
	def cluster = "-k $clusters"
	def ploi = "--ploidy $ploidy"
	def minalt = "--min_alt ${min_alt}"
	def minref = "--min_ref ${min_ref}"
	def maxloci = "--max_loci ${max_loci}"

	def restart = restarts != 'False' ? "--restarts $restarts" : ''
    def commonvariant = common_variants != 'no_common_variants' ? "--common_variants ${common_variants}" : ''
    def with_commonvariant = common_variants != 'no_common_variants' ? "with_common_variant" : 'no_common_variants'

    def knowngenotype = known_genotypes != 'no_known_genotypes' ? "--known_genotypes ${known_genotypes}" : ''
    def with_knowngenotype = known_genotypes != 'no_known_genotypes' ? "with_known_genotype" : 'no_known_genotypes'

    def knowngenotypes_sample = known_genotypes_sample_names != 'False' ? "--known_genotypes_sample_names ${known_genotypes_sample_names}" : ''
	def with_knowngenotypes_sample = known_genotypes_sample_names != 'False' ? "with_knowngenotypes_sample_names" : 'no_knowngenotypes_sample_names'

	def skipremap = skip_remap != 'False' ? "--skip_remap True" : ''
	def ign = ignore != 'False' ? "--ignore True" : ''
	
	def out = "soup+${threads}+${clusters}+${ploidy}+${min_alt}+${min_ref}+${max_loci}+$restarts+${with_commonvariant}+${with_knowngenotype}+${with_knowngenotypes_sample}+${skip_remap}+${ignore}"
    
    """
    singularity exec -B $PWD/souporcell/:/souporcell $projectDir/souporcell_latest.sif souporcell_pipeline.py $bamfile $barcode $fastafile $thread $cluster $ploi $minalt $minref $maxloci $restart $commonvariant $knowngenotype $knowngenotypes_sample $skipremap $ign -o $out
    pwd
    cp -R /root/$out . 
    """  
}


process vireo{
    publishDir "$params.outdir/vireo", mode: 'copy'
    echo true 
    /*
    container "vireo"
    */

    input:
	file celldata
	each ndonor
        
	each vartrixData
	each donorfile
	each genoTag
 
	each noDoublet         
	each nInit 
	each extraDonor
	each extraDonorMode
	each forceLearnGT

	each ASEmode           
	each noPlot            
	each randSeed
	each cellRange
	each callAmbientRNAs
	each nproc


    output:
	path "vireo*"
	

    script:
	def out = "vireo+${ndonor}+${vartrixData}+${donorfile}+${genoTag}+${noDoublet}+${nInit}+${extraDonor}+${extraDonorMode}+${forceLearnGT}+${ASEmode}+${noPlot}+${randSeed}+${cellRange}+${callAmbientRNAs}+${nproc}"

	def cell_data = "-c $celldata"
	def n_donor = "-N $ndonor"
	def vartrix_data = vartrixData != 'no_vartrixData' ? "--vartrixData $vartrixData" : ''
    def donor = donorfile != 'no_donorfile' ? "-d $donorfile" : ''
	def geno_tag = "--genoTag $genoTag"
	def no_doublet = noDoublet != 'False' ? "--noDoublet" : ''
	def n_init = "--nInit $nInit"
    def extra_donor = "--extraDonor $extraDonor"
	def extradonor_mode = extraDonorMode != 'distance' ? "--extraDonorMode $extraDonorMode" : ''
	def learnGT = forceLearnGT != 'False' ? "--forceLearnGT" : ''
	def ase_mode = ASEmode != 'False' ? "--ASEmode" : ''
	def no_plot = noPlot != 'False' ? "--noPlot" : ''
	def random_seed = randSeed != 'none'? "--randSeed $randSeed" : ''
	def cell_range = cellRange != 'all'? "--cellRange $cellRange" : ''
	def call_ambient_rna = callAmbientRNAs != 'False' ? "--callAmbientRNAs" : ''
	def n_proc = "--nproc $nproc"
    
    

    """
    vireo ${cell_data} ${n_donor} ${vartrix_data} $donor ${geno_tag} ${no_doublet} ${n_init} ${extra_donor} ${extradonor_mode} $learnGT ${ase_mode} ${no_plot} ${random_seed} ${cell_range} ${call_ambient_rna} ${n_proc} -o $out  
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


process compare_parameter{
    
	publishDir "$params.outdir/compare", mode: 'copy'
    echo true

    input:
	file demultiplex_result
        
    output:
	file '*.tsv'
	/*
    file "barplot.png"
	file "barplot by group.png"
    */

    script:
    def demux_files = ""
    for(r : demultiplex_result) {
    	demux_files = r + ":" + demux_files
    }

    """
    compareParam.R --file ${demux_files}  
    """
}

workflow demultiplex{
    if (params.demux_tool == "demuxlet"){
	    input_sam = Channel.fromPath(params.sam)
	    input_tag_group = split_input(params.tag_group)
        input_tag_UMI = split_input(params.tag_UMI)
	    input_plp = split_input(params.plp)

	    input_vcf = Channel.fromPath(params.vcf)
	    input_field = split_input(params.field)

	    input_geno_error_offset = split_input(params.geno_error_offset)
        input_geno_error_coeff = split_input(params.geno_error_coeff)
	    input_r2_info= split_input(params.r2_info)
	    input_min_mac = split_input(params.min_mac)
	    input_min_callrate = split_input(params.min_callrate)
	    input_sm = split_input(params.sm)
	    input_sm_list = split_input(params.sm_list)
		input_alpha_list = Channel.value(params.alpha_list)
		input_alpha = split_input(params.alpha)
			
	    input_doublet_prior = split_input(params.doublet_prior)
	    input_sam_verbose = split_input(params.sam_verbose)
		input_vcf_verbose = split_input(params.vcf_verbose)
		
	    input_cap_BQ = split_input(params.cap_BQ)
 	    input_min_BQ = split_input(params.min_BQ)
	    input_min_MQ = split_input(params.min_MQ)
	    input_min_TD = split_input(params.min_TD)
	    input_excl_flag = split_input(params.excl_flag)

		input_group_list = split_input(params.group_list)
	    input_min_total = split_input(params.min_total)
	    input_min_umi = split_input(params.min_umi)    
	    input_min_snp = split_input(params.min_snp)
       
	    demuxlet(input_sam, input_tag_group, input_tag_UMI, input_plp, input_vcf, input_field, input_geno_error_offset, input_geno_error_coeff, input_r2_info, input_min_mac, input_min_callrate, input_sm, input_sm_list, input_alpha_list, input_alpha, input_doublet_prior, input_sam_verbose, input_vcf_verbose, input_cap_BQ, input_min_BQ, input_min_MQ, input_min_TD, input_excl_flag, input_group_list, input_min_total, input_min_umi, input_min_snp)
	    compare_parameter(demuxlet.out.collect())
     
    }
  
    else if (params.demux_tool == "souporcell"){
	    input_bam = channel.value(params.bam)
	    input_barcodes = channel.value(params.barcodes)
	    input_fasta = channel.value(params.fasta)
	    input_threads = split_input(params.threads)
	    input_clusters = split_input(params.clusters)

	    input_ploidy = split_input(params.ploidy)
	    input_min_alt = split_input(params.min_alt)
	    input_min_ref = split_input(params.min_ref)
	    input_max_loci = split_input(params.max_loci)
	    input_restarts = split_input(params.restarts)
  
		input_common_variants = split_input(params.common_variants)
	    input_known_genotypes = split_input(params.known_genotypes)
	    input_known_genotypes_sample_names = split_input(params.known_genotypes_sample_names)
	    input_skip_remap = split_input(params.skip_remap)
	    input_ignore = split_input(params.ignore)
 	   
	    souporcell(input_bam, input_barcodes, input_fasta, input_threads, input_clusters, input_ploidy, input_min_alt, input_min_ref, input_max_loci, input_restarts, input_common_variants, input_known_genotypes, input_known_genotypes_sample_names, input_skip_remap, input_ignore) 	 
		compare_parameter(souporcell.out.collect())
    }
    else if (params.demux_tool == "vireo"){
	    input_celldata = Channel.fromPath(params.celldata)
	    input_ndonor = split_input(params.ndonor)

	    input_vartrixData = split_input(params.vartrixData)
	    input_donorfile = split_input(params.donorfile)
	    input_genoTag = split_input(params.genoTag)

	    input_noDoublet = split_input(params.noDoublet)
	    input_nInit = split_input(params.nInit)
	    input_extraDonor = split_input(params.extraDonor)
	    input_extraDonorMode = split_input(params.extraDonorMode)
	    input_forceLearnGT = split_input(params.forceLearnGT)

		input_ASEmode = split_input(params.ASEmode)
	    input_noPlot = split_input(params.noPlot)
 	    input_randSeed = split_input(params.randSeed)
 	    input_cellRange = split_input(params.cellRange)
 	    input_callAmbientRNAs = split_input(params.callAmbientRNAs)
	    input_nproc = split_input(params.nproc)
     
	    vireo(input_celldata, input_ndonor, input_vartrixData, input_donorfile, input_genoTag, input_noDoublet, input_nInit, input_extraDonor, input_extraDonorMode, input_forceLearnGT, input_ASEmode, input_noPlot, input_randSeed, input_cellRange, input_callAmbientRNAs, input_nproc)
	    compare_parameter(vireo.out.collect())

   }

}

workflow{
	demultiplex()
       
}




