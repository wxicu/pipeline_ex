#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*  
    Note:
	--alpha
	--alpha_list <grid search> becasue --alpha 0 --alpha 0.5 is different from --alpha 0.5 TUT 
    
    echo ${demux_files} > print.txt
    def max_complex_gap = max_complex_gap_freebayes != 'False' ? "--max-complex-gap ${max_complex_gap_freebayes}": ''
    singularity exec -B $PWD/souporcell/:/souporcell $projectDir/souporcell_latest.sif ls $out 
    singularity exec -B $PWD/souporcell/:/souporcell $projectDir/souporcell_latest.sif cp -R $out .
	data = Channel.fromPath(data_path, type: 'dir')	 
    file "barplot.png" file "barplot by group.png"

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
		popscle demuxlet $samfile $taggroup $tagUMI $plpfile $vcffile $fieldinfo ${genoerror_off} ${genoerror_cof} $r2info $minmac $mincallrate $smlist ${sm_list_file} --alpha 0 --alpha $alpha $doubletprior $samverbose $vcfverbose $capBQ $minBQ $minMQ $minTD $exclflag $grouplist $mintotal $minumi $minsnp --out $out
		
		"""  

    else{
   out = "demux+${tag_group}+${tagUMI}+${plp}+${field}+${geno_error_offset}+${geno_error_coeff}+${r2_info}+${min_mac}+${min_callrate}+${sm}+${sm_list}+${alpha_list}+${doublet_prior}+${sam_verbose}+${vcf_verbose}+${cap_BQ}+${min_BQ}+${min_MQ}+${min_TD}+${excl_flag}+${group_list}+${min_total}+${min_umi}+${min_snp}"	
		"""
		popscle demuxlet $sam $taggroup $tagUMI $plpfile $vcffile $fieldinfo ${genoerror_off} ${genoerror_cof} $r2info $minmac $mincallrate $smlist ${sm_list_file} --alpha ${alpha_value} $doubletprior $samverbose $vcfverbose $capBQ $minBQ $minMQ $minTD $exclflag $grouplist $mintotal $minumi $minsnp --out $out
		"""  
	}
	
}

workflow demultiplex_demux{
	main:
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
		
	emit:
	    demuxlet.out.collect()
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

workflow demultiplex_soup{
	main:
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
			
	emit:
	    souporcell.out.collect()
}

process vireo{
    publishDir "$params.outdir/vireo", mode: 'copy'
    echo true 

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

workflow demultiplex_vireo{
	main:
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
	
	emit:
	    vireo.out.collect()
}



process scSplit{
    publishDir "$params.outdir/scSplit", mode: 'copy'
    echo true 

    input:
	each vcf
	each bam
	each barcode
	each tag
	each com
	each ref
	each alt
	each num
	each sub
	each ems
	each dbl
	each vcf_known
        path scSplit_loc
        each sample_geno
   	
    output:
	path "scSplit*"
   
    script:
    def with_commmon = com != 'no_commonData' ? "with_common_snp" : "without_common_snp"
    def with_known = vcf_known != 'no_vcfKnownData' ? "with_known_snp" : "without_known_snp"
    def out = "-o scSplit+${with_commmon}+${num}+${sub}+${ems}+${dbl}+${with_known}"
    def outdir = "scSplit+${with_commmon}+${num}+${sub}+${ems}+${dbl}+${with_known}"
    def vcf_data = "-v $vcf"
    def bam_data = "-i $bam"
    def barcode_data = "-b $barcode"
    def tag_data = "--tag $tag"
    def common_data = com != 'no_commonData' ? "--com $com" : ''
    def num_data = "-n $num"
    def sub_data = "--sub $sub"
    def ems_data = "--ems $ems"
    def dbl_data = dbl != 'False' ? "--dbl $dbl" : ''
    def vcf_known_data = vcf_known != 'no_vcfKnownData' ? "--vcf ${vcf_known}" : ''
    
    if (sample_geno != 'False'){
    """
    mkdir $outdir
    python3 ${scSplit_loc}/scSplit count ${vcf_data} ${bam_data} ${barcode_data} ${common_data} -r $ref -a $alt $out
    python3 ${scSplit_loc}/scSplit run -r ${outdir}/$ref -a ${outdir}/$alt $out ${num_data} ${sub_data} ${ems_data} ${dbl_data} ${vcf_known_data}
    python3 ${scSplit_loc}/scSplit genotype -r ${outdir}/$ref -a ${outdir}/$alt -p ${outdir}/scSplit_P_s_c.csv $out
    """
    }
    else{
    """
    mkdir $outdir
    python3 ${scSplit_loc}/scSplit count ${vcf_data} ${bam_data} ${barcode_data} ${common_data} -r $ref -a $alt $out
    python3 ${scSplit_loc}/scSplit run -r ${outdir}/$ref -a ${outdir}/$alt $out ${num_data} ${sub_data} ${ems_data} ${dbl_data} ${vcf_known_data}
    """
}
    
}

workflow demultiplex_scSplit{
     main:
		input_vcf_scsplit = Channel.fromPath(params.vcfscSplit)
		input_bam_scsplit = split_input(params.bamscSplit)
		input_bar_scsplit = split_input(params.barscSplit)
		input_tag_scsplit = split_input(params.tagscSplit)
		input_com_scsplit = split_input(params.comscSplit)
		input_ref_scsplit = split_input(params.refscSplit)
		input_alt_scsplit = split_input(params.altscSplit)
		input_num_scsplit = split_input(params.numscSplit)
		input_sub_scsplit = split_input(params.subscSplit)
		input_ems_scsplit = split_input(params.emsscSplit)
		input_dbl_scsplit = split_input(params.dblscSplit)
		input_vcf_known_scsplit = split_input(params.vcf_known_scSplit)
                input_scSplit_loc = Channel.from(params.scSplit_loc)
                input_sample_geno = Channel.from(params.sample_geno)
		scSplit(input_vcf_scsplit, input_bam_scsplit, input_bar_scsplit, input_tag_scsplit, input_com_scsplit, input_ref_scsplit, input_alt_scsplit, input_num_scsplit, input_sub_scsplit, input_ems_scsplit, input_dbl_scsplit, input_vcf_known_scsplit,input_scSplit_loc,input_sample_geno)
      emit:
	    scSplit.out.collect()
	


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
    val demultiplex_result
    val tool
    
    output:
    file '*.tsv'

    script:
    def demux_files = ""
    if (tool != "none"){
	for(r : demultiplex_result) {
	    	demux_files = r + ":" + demux_files
    	}
    }
    else{
	demux_files = demultiplex_result
    }
   

    """
    compareParam.R --tool $tool --file ${demux_files} 
    """
}


workflow demultiplex{
    input_tool = Channel.from(params.demux_tool)
    input_compare_file = Channel.from(params.demuxlet_outdir)
	if (params.demux_tool == "demuxlet"){
		demultiplex_demux()
		compare_parameter(demultiplex_demux.out, input_tool)
	}

	else if (params.demux_tool == "souporcell"){
		demultiplex_soup()
		compare_parameter(demultiplex_soup.out, input_tool)
	}

	else if (params.demux_tool == "vireo"){
		demultiplex_vireo()
		compare_parameter(demultiplex_vireo.out, input_tool)
	}
	else if (params.demux_tool == "scSplit"){
		demultiplex_scSplit()
                compare_parameter(demultiplex_scSplit.out, input_tool)
    	}
    	else if (params.demux_tool == "none"){
   		compare_parameter(input_compare_file, input_tool)
    	}
}


process freebayes{
    publishDir "$params.outdir/freebayes", mode: 'copy'
    echo true

    input:
	file bam_freebayes
	file vcf_freebayes
    file ref_freebayes
    val stdin_freebayes
    file targets_freebayes
    file samples_freebayes
    file cnv_map_freebayes
    file trace_freebayes
    file failed_alleles_freebayes
    file variant_input_freebayes
    val only_use_input_alleles_freebayes
    val pvar_freebayes
    val region_freebayes
    val show_reference_repeats_freebayes
    val theta_freebayes
    val ploidy_freebayes
    val pooled_freebayes
    val use_reference_allele_freebayes
    val reference_quality_freebayes
    val no_snps_freebayes
    val no_indels_freebayes
    val no_mnps_freebayes
    val no_complex_freebayes
    val use_best_n_alleles_freebayes
    val use_duplicate_reads_freebayes
    val min_mapping_quality_freebayes
    val min_base_quality_freebayes
    val mismatch_base_quality_threshold_freebayes
    val read_mismatch_limit_freebayes
    val left_align_indels_freebayes
    val read_max_mismatch_fraction_freebayes
    val read_snp_limit_freebayes
    val read_indel_limit_freebayes
    val no_filters_freebayes
    val indel_exclusion_window_freebayes
    val min_alternate_fraction_freebayes
    val min_alternate_count_freebayes
    val min_alternate_qsum_freebayes
    val min_alternate_total_freebayes
    val min_coverage_freebayes
    val no_ewens_priors_freebayes
    val no_population_priors_freebayes
    val hwe_priors_freebayes
    val binomial_obs_priors_freebayes
    val allele_balance_priors_freebayes
    val site_selection_max_iterations_freebayes
    val genotyping_max_iterations_freebayes
    val posterior_integration_limits_freebayes
    val no_permute_freebayes
    val exclude_unobserved_genotypes_freebayes
    val genotype_variant_threshold_freebayes
    val use_mapping_quality_freebayes
    val read_dependence_factor_freebayes
    val no_marginals_freebayes
    val debug_freebayes
    val dd_freebayes
    val gvcf_freebayes
    val gvcf_chunk_freebayes
    val gvcf_dont_use_chunk_freebayes
    file haplotype_basis_alleles_freebayes
    val report_all_haplotype_alleles_freebayes
    val report_monomorphic_freebayes
    val strict_vcf_freebayes
    val pooled_discrete_freebayes
    val pooled_continuous_freebayes 
    val haplotype_length_freebayes 
    val min_repeat_size_freebayes
    val min_repeat_entropy_freebayes 
    val no_partial_observations_freebayes 
    val min_supporting_allele_qsum_freebayes 
    val min_supporting_mapping_qsum_freebayes
    val standard_filters_freebayes
    val limit_coverage_freebayes
    val skip_coverage_freebayes 
    file observation_bias_freebayes
    val base_quality_cap_freebayes  
    val prob_contamination_freebayes
    val legacy_gls_freebayes 
    file contamination_estimates_freebayes 
    val report_genotype_likelihood_max_freebayes
    val genotyping_max_banddepth_freebayes 
    val harmonic_indel_quality_freebayes 
    val genotype_qualities_freebayes
    val throw_away_snp_obs_freebayes
    val throw_away_indel_obs_freebayes
    val throw_away_mnps_obs_freebayes
    val throw_away_complex_obs_freebayes


    output:
	file '*.vcf'

    script:
    def stdin = stdin_freebayes != 'False' ? "--stdin" : ''
    def targets = targets_freebayes.name != 'no_targets' ? "--targets ${targets_freebayes}" : ''
    def samples = samples_freebayes.name != 'no_samples' ? "--samples ${samples_freebayes}" : ''
    def cnv_map = cnv_map_freebayes.name != 'no_cnv_map' ? "--cnv-map ${cnv_map_freebayes}" : ''
    def trace = trace_freebayes.name != 'no_trace' ? "--trace ${trace_freebayes_freebayes}" : ''
    def failed_alleles = failed_alleles_freebayes.name != 'no_failed_alleles' ? "--failed-alleles ${failed_alleles_freebayes}" : ''
    def variant_input = variant_input_freebayes.name != 'no_variant_input' ? "--variant-input ${variant_input_freebayes}" : ''
    def only_use_input_alleles = only_use_input_alleles_freebayes != 'False' ? "--only-use-input-alleles" : ''
    def pvar = "--pvar ${pvar_freebayes}" 
    def region = region_freebayes != 'False' ? "--region ${region_freebayes}" : ''
    def show_reference_repeats = show_reference_repeats_freebayes != 'False' ? "--show-reference-repeats" : ''
    def theta = "--theta ${theta_freebayes}"
    def ploidy = "--ploidy ${ploidy_freebayes}" 
    def pool = pooled_freebayes != 'False' ? "--pooled" : ''
    def use_reference_allele = use_reference_allele_freebayes != 'False' ? "--use-reference-allele" : ''
    def reference_quality = "--reference-quality ${reference_quality_freebayes}"  
    def no_snps = no_snps_freebayes != 'False' ? "--no-snps" : ''
    def no_indels = no_indels_freebayes != 'False' ? "--no-indels" : ''
    def no_mnps = no_mnps_freebayes != 'False' ? "--no-mnps" : ''
    def no_complex = no_complex_freebayes != 'False' ? "--no-complex" : ''
    def use_best_n_alleles = "--use-best-n-alleles ${use_best_n_alleles_freebayes}"

    def use_duplicate_reads = use_duplicate_reads_freebayes != 'false' ? "--use-duplicate-reads" : ''
    def min_mapping_quality = "--min-mapping-quality ${min_mapping_quality_freebayes}"
    def min_base_quality = "--min-base-quality ${min_base_quality_freebayes}" 
    def mismatch_base_quality_threshold = "--mismatch-base-quality-threshold ${mismatch_base_quality_threshold_freebayes}"
    def read_mismatch_limit = read_mismatch_limit_freebayes != 'False' ? "-read-mismatch-limit ${read_mismatch_limit_freebayes}":''
    def left_align_indels = left_align_indels_freebayes != 'false' ? "--dont-left-align-indels" : ''
    def read_max_mismatch_fraction = "--read-max-mismatch-fraction ${read_max_mismatch_fraction_freebayes}"
    def read_snp_limit = read_snp_limit_freebayes != 'False' ? "--read-snp-limit ${read_snp_limit_freebayes}":''
    def read_indel_limit = read_indel_limit_freebayes != 'False' ? "--read-indel-limit ${read_indel_limit_freebayes}" :''
    def no_filters = no_filters_freebayes != 'False' ? "--no-filters" : ''
    def indel_exclusion_window = "--indel-exclusion-window ${indel_exclusion_window_freebayes}"
    def min_alternate_fraction = "--min-alternate-fraction ${min_alternate_fraction_freebayes}"
    def min_alternate_count = "--min-alternate-count ${min_alternate_count_freebayes}"
    def min_alternate_qsum = "--min-alternate-qsum ${min_alternate_qsum_freebayes}"
    def min_alternate_total = "--min-alternate-total ${min_alternate_total_freebayes}"
    def min_coverage = "--min-coverage ${min_coverage_freebayes}"
    def no_ewens_priors = no_ewens_priors_freebayes != 'False' ? "--no-ewens-priors" : ''
    def no_population_priors = no_population_priors_freebayes != 'False' ? "--no-population-priors" : ''
    def hwe_priors = hwe_priors_freebayes != 'False' ? "--hwe-priors-off" : ''
    def binomial_obs_priors = binomial_obs_priors_freebayes != 'False' ? "--binomial-obs-priors-off" : ''
    def allele_balance_priors = allele_balance_priors_freebayes != 'False' ? "--allele-balance-priors-off" : ''
    def site_selection_max_iterations = "--site-selection-max-iterations ${site_selection_max_iterations_freebayes}"
    def genotyping_max_iterations = "--genotyping-max-iterations ${genotyping_max_iterations_freebayes}" 
    def posterior_integration_limits = "--posterior-integration-limits ${posterior_integration_limits_freebayes}" 
    def no_permute = no_permute_freebayes != 'False' ? "--no-permute" : ''
    def exclude_unobserved_genotypes = exclude_unobserved_genotypes_freebayes != 'False' ? "--exclude-unobserved-genotypes" : ''
    def genotype_variant_threshold = genotype_variant_threshold_freebayes != 'False' ? "--genotype-variant-threshold ${genotype_variant_threshold_freebayes}":''
    def use_mapping_quality = use_mapping_quality_freebayes != 'False' ? "--use-mapping-quality" : ''
    def read_dependence_factor = "--read-dependence-factor ${read_dependence_factor_freebayes}"
    def no_marginals = no_marginals_freebayes != 'False' ? "--no-marginals" : ''
    def debug = debug_freebayes != 'False' ? "--debug" : ''
    def dd = dd_freebayes != 'False' ? "-dd" : ''
    def gvcf = gvcf_freebayes !='False' ? "--gvcf":''
    def gvcf_chunk = gvcf_chunk_freebayes !='False' ? "--gvcf-chunk ${gvcf_chunk_freebayes}":''
    def gvcf_dont_use_chunk = gvcf_dont_use_chunk_freebayes !='false' ?"--gvcf-dont-use-chunk ${gvcf_dont_use_chunk_freebayes}":''
    def haplotype_basis_alleles = haplotype_basis_alleles_freebayes.name != 'no_haplotype_basis_alleles' ? "--haplotype-basis-alleles ${haplotype_basis_alleles_freebayes}":''
    def report_all_haplotype_alleles = report_all_haplotype_alleles_freebayes !='False' ? "--report-all-haplotype-alleles":''
    def report_monomorphic = report_monomorphic_freebayes !='False' ? "--report-monomorphic":''
    def strict_vcf = strict_vcf_freebayes!='False' ? "--strict-vcf":''
    def pooled_discrete = pooled_discrete_freebayes!='False' ? "--pooled-discrete":''
    def pooled_continuous = pooled_continuous_freebayes !='False' ? "--pooled-continuous":''
    def haplotype_length = haplotype_length_freebayes = "--haplotype-length ${haplotype_length_freebayes}"
    def min_repeat_size = "--min-repeat-size ${min_repeat_size_freebayes}"
    def min_repeat_entropy = "--min-repeat-entropy ${min_repeat_entropy_freebayes}"
    def no_partial_observations = no_partial_observations_freebayes !='False' ? "--no-partial-observations":''
    def min_supporting_allele_qsum = "--min-supporting-allele-qsum ${min_supporting_allele_qsum_freebayes}"
    def min_supporting_mapping_qsum = "--min-supporting-mapping-qsum ${min_supporting_mapping_qsum_freebayes}"
    def standard_filters = standard_filters_freebayes!='False' ? "--standard-filters":''
    def limit_coverage = limit_coverage_freebayes!='False' ? "--limit-coverage ${limit_coverage_freebayes}":''
    def skip_coverage = skip_coverage_freebayes !='False' ? "--skip-coverage ${skip_coverage_freebayes}":''
    def observation_bias = observation_bias_freebayes.name != 'no_observation_bias' ? "--observation-bias ${observation_bias_freebayes}":''
    def base_quality_cap = base_quality_cap_freebayes !='False' ? "--base-quality-cap ${base_quality_cap_freebayes}":''
    def prob_contamination = prob_contamination_freebayes = "--prob-contamination ${prob_contamination_freebayes}"
    def legacy_gls = legacy_gls_freebayes !='False' ? "--legacy-gls" :''
    def contamination_estimates = contamination_estimates_freebayes.name != 'no_contamination_estimates' ? "--contamination-estimates ${contamination_estimates_freebayes}":''
    def report_genotype_likelihood_max = report_genotype_likelihood_max_freebayes !='False' ? "--report-genotype-likelihood-max":''
    def genotyping_max_banddepth = "--genotyping-max-banddepth ${genotyping_max_banddepth_freebayes}"
    def harmonic_indel_quality = harmonic_indel_quality_freebayes !='False' ? "--harmonic-indel-quality" :''
    def genotype_qualities = genotype_qualities_freebayes!='False' ? "--genotype-qualities" :''
    def throw_away_snp_obs = throw_away_snp_obs_freebayes!='False' ? "--throw-away-snp-obs " :''
    def throw_away_indel_obs = throw_away_indel_obs_freebayes!='False' ? "--throw-away-indel-obs " :''
    def throw_away_mnps_obs = throw_away_mnps_obs_freebayes!='False' ? "--throw-away-mnps-obs " :''
    def throw_away_complex_obs = throw_away_complex_obs_freebayes!='False' ? "--throw-away-complex-obs " :''

    """
    freebayes -f ${ref_freebayes} ${bam_freebayes} $stdin $targets $samples ${cnv_map} $trace ${failed_alleles} ${variant_input} ${only_use_input_alleles} $pvar $region ${show_reference_repeats} $theta $ploidy $pool ${use_reference_allele} ${reference_quality} ${no_snps} ${no_indels} ${no_mnps} ${no_complex} ${use_best_n_alleles} ${use_duplicate_reads} ${min_mapping_quality} ${min_base_quality} ${mismatch_base_quality_threshold} ${read_mismatch_limit} ${left_align_indels} ${read_max_mismatch_fraction} $read_snp_limit ${read_indel_limit} ${no_filters} ${indel_exclusion_window} ${min_alternate_fraction} ${min_alternate_count} ${min_alternate_qsum} ${min_alternate_total} ${min_coverage} ${no_ewens_priors} ${no_population_priors} ${hwe_priors} ${binomial_obs_priors} ${allele_balance_priors} ${site_selection_max_iterations} ${genotyping_max_iterations} ${posterior_integration_limits} ${no_permute} ${exclude_unobserved_genotypes} ${genotype_variant_threshold} ${use_mapping_quality} ${read_dependence_factor} ${no_marginals} $debug $dd ${gvcf} ${gvcf_chunk} ${gvcf_dont_use_chunk} ${haplotype_basis_alleles} ${report_all_haplotype_alleles} ${report_monomorphic} ${strict_vcf} ${pooled_discrete} ${pooled_continuous} ${haplotype_length} ${min_repeat_size} ${min_repeat_entropy} ${no_partial_observations} ${min_supporting_allele_qsum} ${min_supporting_mapping_qsum} ${standard_filters} ${limit_coverage} ${skip_coverage} ${observation_bias} ${base_quality_cap} ${prob_contamination} ${legacy_gls} ${contamination_estimates} ${report_genotype_likelihood_max} ${genotyping_max_banddepth} ${harmonic_indel_quality} ${genotype_qualities} ${throw_away_snp_obs} ${throw_away_indel_obs} ${throw_away_mnps_obs} ${throw_away_complex_obs} > ${vcf_freebayes} 
    """
}


workflow variant_freebayes{
	input_bam_freebayes = Channel.fromPath(params.bam_freebayes)
	input_vcf_freebayes = channel.fromPath(params.vcf_freebayes)
	input_fasta_reference = Channel.fromPath(params.fasta_reference)
	input_stdin_freebayes = Channel.value(params.stdin)
	input_targets_freebayes = Channel.fromPath(params.targets)
	input_samples_freebayes = Channel.fromPath(params.samples)
	input_cnv_map_freebayes = Channel.fromPath(params.cnv_map)
	input_trace_freebayes = Channel.fromPath(params.trace)
	input_failed_alleles_freebayes = Channel.fromPath(params.failed_alleles)
	input_variant_input_freebayes = Channel.fromPath(params.variant_input)
    input_only_use_input_alleles_freebayes = Channel.value(params.only_use_input_alleles)
	input_pvar_freebayes = channel.value(params.pvar)
    input_region_freebayes = channel.value(params.region)
	input_show_reference_repeats_freebayes = channel.value(params.show_reference_repeats)
    input_theta_freebayes = channel.value(params.theta)
	input_ploidy_freebayes = channel.value(params.ploidy)
	input_pooled_freebayes = channel.value(params.pooled)
	input_use_reference_allele_freebayes = channel.value(params.use_reference_allele)
	input_reference_quality_freebayes = channel.value(params.reference_quality)
	input_no_snps_freebayes = channel.value(params.no_snps)
    input_no_indels_freebayes = channel.value(params.no_indels)
    input_no_mnps_freebayes = channel.value(params.no_mnps)
    input_no_complex_freebayes = channel.value(params.no_complex)
    input_use_best_n_alleles_freebayes = channel.value(params.use_best_n_alleles)
    input_use_duplicate_reads_freebayes = channel.value(params.use_duplicate_reads)
    input_min_mapping_quality_freebayes = channel.value(params.min_mapping_quality)
    input_min_base_quality_freebayes = channel.value(params.min_base_quality)
    input_mismatch_base_quality_threshold_freebayes = channel.value(params.mismatch_base_quality_threshold)
    input_read_mismatch_limit_freebayes = channel.value(params.read_mismatch_limit)    
    input_left_align_indels_freebayes = channel.value(params.left_align_indels)
    input_read_max_mismatch_fraction_freebayes = channel.value(params.read_max_mismatch_fraction)
    input_read_snp_limit_freebayes = channel.value(params.read_snp_limit)
    input_read_indel_limit_freebayes = channel.value(params.read_indel_limit)
    input_no_filters_freebayes = channel.value(params.no_filters)
    input_indel_exclusion_window_freebayes = channel.value(params.indel_exclusion_window)
    input_min_alternate_fraction_freebayes = channel.value(params.min_alternate_fraction)
    input_min_alternate_count_freebayes = channel.value(params.min_alternate_count)
    input_min_alternate_qsum_freebayes = channel.value(params.min_alternate_qsum)
    input_min_alternate_total_freebayes = channel.value(params.min_alternate_total)
    input_min_coverage_freebayes = channel.value(params.min_coverage)
    input_no_ewens_priors_freebayes = channel.value(params.no_ewens_priors)
    input_no_population_priors_freebayes = channel.value(params.no_population_priors)
    input_hwe_priors_freebayes = channel.value(params.hwe_priors)
    input_binomial_obs_priors_freebayes = channel.value(params.binomial_obs_priors)
    input_allele_balance_priors_freebayes = channel.value(params.allele_balance_priors)
    input_site_selection_max_iterations_freebayes = channel.value(params.site_selection_max_iterations)
    input_genotyping_max_iterations_freebayes = channel.value(params.genotyping_max_iterations)
    input_posterior_integration_limits_freebayes = channel.value(params.posterior_integration_limits)
    input_no_permute_freebayes = channel.value(params.no_permute)
    input_exclude_unobserved_genotypes_freebayes = channel.value(params.exclude_unobserved_genotypes)
    input_genotype_variant_threshold_freebayes = channel.value(params.genotype_variant_threshold)
    input_use_mapping_quality_freebayes = channel.value(params.use_mapping_quality)
    input_read_dependence_factor_freebayes = channel.value(params.read_dependence_factor)
    input_no_marginals_freebayes = channel.value(params.no_marginals)
    input_debug_freebayes = channel.value(params.debug)
    input_dd_freebayes = channel.value(params.dd)
    input_gvcf_freebayes = channel.value(params.gvcf)
    input_gvcf_chunk_freebayes = channel.value(params.gvcf_chunk)
    input_gvcf_dont_use_chunk_freebayes = channel.value(params.gvcf_dont_use_chunk)
    input_haplotype_basis_alleles_freebayes = channel.fromPath(params.haplotype_basis_alleles)
    input_report_all_haplotype_alleles_freebayes = channel.value(params.report_all_haplotype_alleles)
    input_report_monomorphic_freebayes = channel.value(params.report_monomorphic)
    input_strict_vcf_freebayes = channel.value(params.strict_vcf)
    input_pooled_discrete_freebayes = channel.value(params.pooled_discrete)
    input_pooled_continuous_freebayes = channel.value(params.pooled_continuous)
    input_haplotype_length_freebayes = channel.value(params.haplotype_length)
    input_min_repeat_size_freebayes = channel.value(params.min_repeat_size)
    input_min_repeat_entropy_freebayes = channel.value(params.min_repeat_entropy)
    input_no_partial_observations_freebayes = channel.value(params.no_partial_observations)
    input_min_supporting_allele_qsum_freebayes = channel.value(params.min_supporting_allele_qsum)
    input_min_supporting_mapping_qsum_freebayes = channel.value(params.min_supporting_mapping_qsum)
    input_standard_filters_freebayes = channel.value(params.standard_filters)
    input_limit_coverage_freebayes = channel.value(params.limit_coverage)
    input_skip_coverage_freebayes = channel.value(params.skip_coverage)
    input_observation_bias_freebayes = channel.fromPath(params.observation_bias)
    input_base_quality_cap_freebayes = channel.value(params.base_quality_cap)
    input_prob_contamination_freebayes = channel.value(params.prob_contamination)
    input_legacy_gls_freebayes = channel.value(params.legacy_gls)
    input_contamination_estimates_freebayes = channel.fromPath(params.contamination_estimates)
    input_report_genotype_likelihood_max_freebayes = channel.value(params.report_genotype_likelihood_max)
    input_genotyping_max_banddepth_freebayes = channel.value(params.genotyping_max_banddepth)
    input_harmonic_indel_quality_freebayes = channel.value(params.harmonic_indel_quality)
    input_genotype_qualities_freebayes = channel.value(params.genotype_qualities)
    input_throw_away_snp_obs = channel.value(params.throw_away_snp_obs)
    input_throw_away_indel_obs = channel.value(params.throw_away_indel_obs)
    input_throw_away_mnps_obs = channel.value(params.throw_away_mnps_obs)
    input_throw_away_complex_obs = channel.value(params.throw_away_complex_obs)
 

    freebayes(input_bam_freebayes,input_vcf_freebayes,input_fasta_reference,input_stdin_freebayes,input_targets_freebayes,input_samples_freebayes,input_cnv_map_freebayes,input_trace_freebayes,input_failed_alleles_freebayes,input_variant_input_freebayes,input_only_use_input_alleles_freebayes,input_pvar_freebayes,input_region_freebayes,input_show_reference_repeats_freebayes,input_theta_freebayes,input_ploidy_freebayes,input_pooled_freebayes,input_use_reference_allele_freebayes,input_reference_quality_freebayes,input_no_snps_freebayes,input_no_indels_freebayes,input_no_mnps_freebayes,input_no_complex_freebayes,input_use_best_n_alleles_freebayes,input_use_duplicate_reads_freebayes,input_min_mapping_quality_freebayes,input_min_base_quality_freebayes,input_mismatch_base_quality_threshold_freebayes,input_read_mismatch_limit_freebayes,input_left_align_indels_freebayes,input_read_max_mismatch_fraction_freebayes,input_read_snp_limit_freebayes, input_read_indel_limit_freebayes,input_no_filters_freebayes,input_indel_exclusion_window_freebayes,input_min_alternate_fraction_freebayes,input_min_alternate_count_freebayes,input_min_alternate_qsum_freebayes,input_min_alternate_total_freebayes,input_min_coverage_freebayes,input_no_ewens_priors_freebayes,input_no_population_priors_freebayes,input_hwe_priors_freebayes,input_binomial_obs_priors_freebayes,input_allele_balance_priors_freebayes,input_site_selection_max_iterations_freebayes,input_genotyping_max_iterations_freebayes,input_posterior_integration_limits_freebayes,input_no_permute_freebayes,input_exclude_unobserved_genotypes_freebayes,input_genotype_variant_threshold_freebayes,input_use_mapping_quality_freebayes,input_read_dependence_factor_freebayes,input_no_marginals_freebayes,input_debug_freebayes,input_dd_freebayes,input_gvcf_freebayes,input_gvcf_chunk_freebayes,input_gvcf_dont_use_chunk_freebayes,input_haplotype_basis_alleles_freebayes,input_report_all_haplotype_alleles_freebayes,input_report_monomorphic_freebayes,input_strict_vcf_freebayes,input_pooled_discrete_freebayes,input_pooled_continuous_freebayes,input_haplotype_length_freebayes,input_min_repeat_size_freebayes,input_min_repeat_entropy_freebayes, input_no_partial_observations_freebayes, input_min_supporting_allele_qsum_freebayes, input_min_supporting_mapping_qsum_freebayes, input_standard_filters_freebayes, input_limit_coverage_freebayes, input_skip_coverage_freebayes, input_observation_bias_freebayes, input_base_quality_cap_freebayes, input_prob_contamination_freebayes, input_legacy_gls_freebayes, input_contamination_estimates_freebayes,input_report_genotype_likelihood_max_freebayes,input_genotyping_max_banddepth_freebayes,input_harmonic_indel_quality_freebayes, input_genotype_qualities_freebayes, input_throw_away_snp_obs, input_throw_away_indel_obs, input_throw_away_mnps_obs, input_throw_away_complex_obs)

}
process cellSNP{
    publishDir "$params.outdir/cellSNP", mode: 'copy'
    echo true

    input:
	
    val samFile_cellSNP
    file samFileList_cellSNP
    file regionsVCF_cellSNP
    file targetsVCF_cellSNP
    file barcodeFile_cellSNP
    file sampleList_cellSNP
    val sampleIDs_cellSNP
    val genotype_cellSNP
    val gzip_cellSNP
    val printSkipSNPs_cellSNP
    val nproc_cellSNP
    file refseq_cellSNP
    val chrom_cellSNP
    val cellTAG_cellSNP
    val UMItag_cellSNP
    val minCOUNT_cellSNP
    val minMAF_cellSNP
    val doubletGL_cellSNP
    val inclFLAG_cellSNP
    val exclFLAG_cellSNP
    val minLEN_cellSNP
    val minMAPQ_cellSNP
    val maxDEPTH_cellSNP
    val countORPHAN_cellSNP
    val outDir_cellSNP


    output:
    file "cellSNP*"

    script:
    def samFile = samFile_cellSNP != 'no_samFile' ? "--samFile ${samFile_cellSNP}" : ''
    def samFileList = samFileList_cellSNP.name != 'no_samFileList' ? "--samFileList ${samFileList_cellSNP}" : ''
    def regionsVCF = regionsVCF_cellSNP.name != 'no_regionsVCF' ? "--regionsVCF ${regionsVCF_cellSNP}" : ''
    def targetsVCF = targetsVCF_cellSNP.name != 'no_targetsVCF' ? "--targetsVCF ${targetsVCF_cellSNP}" : ''
    def barcodeFile = barcodeFile_cellSNP.name != 'no_barcodeFile' ? "--barcodeFile ${barcodeFile_cellSNP}" : ''
    def sampleList = sampleList_cellSNP.name != 'no_sampleList' ? "--sampleList ${sampleList_cellSNP}" : ''
    def sampleIDs = sampleIDs_cellSNP != 'no_sampleIDs' ? "--sampleIDs ${sampleIDs_cellSNP}" : ''
    def genotype = genotype_cellSNP != 'False' ? "--genotype" : ''
    def gzip = gzip_cellSNP != 'False' ? "--gzip" : ''
    def printSkipSNPs = printSkipSNPs_cellSNP != 'False' ? "--printSkipSNPs" : ''
    def nproc = nproc_cellSNP != 'False' ? "--nproc ${nproc_cellSNP}" : ''
    def refseq = refseq_cellSNP.name != 'no_refseq' ? "--refseq ${refseq_cellSNP}" : ''
    def chrom = chrom_cellSNP != 'False' ? "--chrom ${chrom_cellSNP}" : ''
    def cellTAG = "--cellTAG ${cellTAG_cellSNP}" 
    def UMItag = "--UMItag ${UMItag_cellSNP}" 
    def minCOUNT = "--minCOUNT ${minCOUNT_cellSNP}" 
    def minMAF = "--minMAF ${minMAF_cellSNP}" 
    def doubletGL = doubletGL_cellSNP != 'False' ? "--doubletGL" : ''
    def inclFLAG = inclFLAG_cellSNP != 'False' ? "--inclFLAG ${inclFLAG_cellSNP}" : ''
    def exclFLAG = exclFLAG_cellSNP != 'False' ? "--exclFLAG ${exclFLAG_cellSNP}" : ''
    def minLEN = "--minLEN ${minLEN_cellSNP}" 
    def minMAPQ = "--minMAPQ ${minMAPQ_cellSNP}" 
    def maxDEPTH = "--maxDEPTH ${maxDEPTH_cellSNP}" 
    def countORPHAN = countORPHAN_cellSNP != 'False' ? "--countORPHAN" : ''
    def outDir = "--outDir ${outDir_cellSNP}" 

    """
    cellsnp-lite $samFile $samFileList $regionsVCF $targetsVCF $barcodeFile $sampleList $sampleIDs $genotype $gzip $printSkipSNPs $nproc $refseq $chrom $cellTAG $UMItag $minCOUNT $minMAF $doubletGL $inclFLAG $exclFLAG $minLEN $minMAPQ $maxDEPTH $countORPHAN $outDir
    """
}


workflow variant_cellSNP{
    input_samFile = channel.value(params.samFile)
    input_samFileList = channel.fromPath(params.samFileList)
    input_regionsVCF = channel.fromPath(params.regionsVCF)
    input_targetsVCF =  channel.fromPath(params.targetsVCF)
    input_barcodeFile = channel.fromPath(params.barcodeFile)
    input_sampleList = channel.fromPath(params.sampleList)
    input_sampleIDs = channel.value(params.sampleIDs)
    input_genotype_cellSNP = channel.value(params.genotype_cellSNP)
    input_gzip_cellSNP = channel.value(params.gzip_cellSNP)
    input_printSkipSNPs = channel.value(params.printSkipSNPs)
    input_nproc_cellSNP = channel.value(params.nproc_cellSNP)
    input_refseq_cellSNP = channel.fromPath(params.refseq_cellSNP)
    input_chrom = channel.value(params.chrom)
    input_cellTAG = channel.value(params.cellTAG)
    input_UMItag = channel.value(params.UMItag)
    input_minCOUNT = channel.value(params.minCOUNT)
    input_minMAF = channel.value(params.minMAF)
    input_doubletGL = channel.value(params.doubletGL)
    input_inclFLAG = channel.value(params.inclFLAG)
    input_exclFLAG = channel.value(params.exclFLAG)
    input_minLEN = channel.value(params.minLEN)
    input_minMAPQ = channel.value(params.minMAPQ)
    input_maxDEPTH = channel.value(params.maxDEPTH)
    input_countORPHAN = channel.value(params.countORPHAN)
    input_outDir_cellSNP = channel.value(params.outDir_cellSNP)

    cellSNP(input_samFile, input_samFileList, input_regionsVCF, input_targetsVCF, input_barcodeFile, input_sampleList, input_sampleIDs, input_genotype_cellSNP, input_gzip_cellSNP, input_printSkipSNPs, input_nproc_cellSNP, input_refseq_cellSNP, input_chrom, input_cellTAG, input_UMItag, input_minCOUNT, input_minMAF, input_doubletGL, input_inclFLAG, input_exclFLAG, input_minLEN, input_minMAPQ, input_maxDEPTH, input_countORPHAN,input_outDir_cellSNP)


}
workflow call_variant{
	if (params.variant_tool == "freebayes"){
		variant_freebayes()
	}
    if (params.variant_tool == "cellSNP"){
		variant_cellSNP()
	}
}

workflow{
	call_variant()  
}

