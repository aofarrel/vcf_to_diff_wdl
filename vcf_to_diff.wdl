version 1.0

task make_mask_and_diff {
	# This is the version currently used by myco!
	#
	# This combines the creation of the bed graph histogram mask file and the
	# creation of the diff file into one WDL task. Sometimes, in WDL, it is
	# easier to combine tasks to avoid shenanigans with the Array[File] type.
	# The masking part of this task is based on github.com/aofarrel/mask-by-coverage
	input {
		File bam
		Boolean force_diff = false
		Boolean histograms = false
		Float max_ratio_low_coverage_sites_per_sample # for sra, default was 0.05
		Int min_coverage_per_site
		File? tbmf
		File vcf

		# runtime attributes
		Int addldisk = 10
		Int cpu      = 8
		Int retries  = 1
		Int memory   = 16
		Int preempt  = 1
	}
	String basename_bam = basename(bam, ".bam")
	String basename_vcf = basename(vcf, ".vcf")
	Int finalDiskSize = ceil(size(bam, "GB")*2) + ceil(size(vcf, "GB")*2) + addldisk

	parameter_meta {
		bam: "BAM file for this sample"
		force_diff: "Output a diff file even if sample is discarded for being below min_proportion_low_coverage_per_sample"
		histograms: "Generate histogram output"
		max_ratio_low_coverage_sites_per_sample: "If over this percent (0.5 = 50%) of a sample's sites get masked due to being below min_coverage_per_site, discard the entire sample"
		min_coverage_per_site: "Positions with coverage below this value will be masked in diff files"
		tbmf: "BED file of regions of the genome you always want to mask (default: R00000039_repregions.bed)"
		vcf: "VCF file for this sample"
	}
	
	command <<<
	set -eux pipefail
	start=$(date +%s)

	# We want the mask file the user input, if it exists, to be the mask file, else
	# fall back on a default mask file that exists in the Docker image already.
	# We cannot use WDL built-in select_first, or else the user inputting a mask file
	# will result in WDL looking for the literal gs:// URI rather than while the file
	# is localized. Different WDL executors localize files to different places, so the
	# following workaround, while goofy, seems to be the most robust.
	if [[ "~{tbmf}" = "" ]]
	then
		mask="/mask/R00000039_repregions.bed"
	else
		mask="~{tbmf}"
	fi
	
	echo "Copying bam..."
	cp ~{bam} .
	
	echo "Sorting bam..."
	samtools sort -u ~{basename_bam}.bam > sorted_u_~{basename_bam}.bam
	
	echo "Calculating coverage..."
	bedtools genomecov -ibam sorted_u_~{basename_bam}.bam -bga | \
		awk '$4 < ~{min_coverage_per_site}' > \
		~{basename_bam}_below_~{min_coverage_per_site}x_coverage.bedgraph
	
	if [[ "~{histograms}" = "true" ]]
	then
		echo "Generating histograms..."
		bedtools genomecov -ibam sorted_u_~{basename_bam}.bam > histogram.txt
	fi
	
	echo "Pulling diff script..."
	wget https://raw.githubusercontent.com/aofarrel/parsevcf/1.3.1/vcf_to_diff_script.py
	echo "Running script..."
	python3 vcf_to_diff_script.py -v ~{vcf} \
	-d . \
	-tbmf ${mask} \
	-bed ~{basename_bam}_below_~{min_coverage_per_site}x_coverage.bedgraph \
	-cd ~{min_coverage_per_site}
	
	# if the sample has too many low coverage sites, throw it out
	this_files_info=$(awk -v file_to_check="~{basename_vcf}.diff" '$1 == file_to_check' "~{basename_vcf}.report")
	if [[ ! "$this_files_info" = "" ]]
	then
		# okay, we have information about this file. is it above the removal threshold?
		echo "$this_files_info" > temp
		amount_low_coverage=$(cut -f2 temp)
		
		# account for very tiny numbers
		if [[ "$amount_low_coverage" == *"e"* ]]
		then
			echo "Scientific notation detected, so it's likely this sample is very much passing."
			echo "PASS" >> ERROR
			exit 0
		fi
		
		percent_low_coverage=$(echo "$amount_low_coverage"*100 | bc)
		echo "$percent_low_coverage percent of ~{basename_vcf} is below ~{min_coverage_per_site}x coverage."
		
		# piping an inequality to `bc` will return 0 if false, 1 if true
		is_bigger=$(echo "$amount_low_coverage>~{max_ratio_low_coverage_sites_per_sample}" | bc)
		if [[ $is_bigger == 0 ]]
		then
			# amount of low coverage is BELOW the removal threshold: sample passes
			echo "PASS" >> ERROR
	
		else
			# amount of low coverage is ABOVE the removal threshold: sample fails
			if [[ "~{force_diff}" == "false" ]]
			then
				rm "~{basename_vcf}.diff"
			fi
			pretty_percent=$(printf "%0.2f" "$percent_low_coverage")
			echo FAILURE - "$pretty_percent""%" is above "~{max_ratio_low_coverage_sites_per_sample}""%" cutoff
			echo VCF2DIFF_"$pretty_percent"_PCT_BELOW_"~{min_coverage_per_site}"x_COVERAGE >> ERROR
		fi
	fi
			
	end=$(date +%s)
	seconds=$(echo "$end - $start" | bc)
	minutes=$(echo "$seconds" / 60 | bc)
	echo "Finished in about $minutes minutes ($seconds sec))"
	ls -lha
	>>>

	runtime {
		cpu: cpu
		docker: "ashedpotatoes/sranwrp:1.1.15"
		disks: "local-disk " + finalDiskSize + " HDD"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	meta {
		author: "Lily Karim (WDLization by Ash O'Farrell)"
	}

	output {
		File mask_file = basename_bam+"_below_"+min_coverage_per_site+"x_coverage.bedgraph"  # !StringCoercion
		File? diff = basename_vcf+".diff"
		File? report = basename_vcf+".report"
		File? histogram = "histogram.txt"
		String errorcode = read_string("ERROR")
	}
}

task make_diff_from_vcf_and_mask {
	input {
		File vcf
		File? tbmf
		Float max_ratio_low_coverage_sites_per_sample = 0.05
		Boolean force_diff = false
		Int min_coverage_per_site = 10
		File bedgraph
		
		# runtime attributes
		Int addldisk = 10
		Int cpu      = 8
		Int retries  = 1
		Int memory   = 16
		Int preempt  = 1
	}
	String basename_vcf = basename(vcf, ".vcf")
	Int finalDiskSize = ceil(size(bedgraph, "GB")*2) + ceil(size(vcf, "GB")*2) + addldisk
	
	command <<<
	set -eux pipefail

	# We want the mask file the user input, if it exists, to be the mask file, else
	# fall back on a default mask file that exists in the Docker image already.
	# We cannot use WDL built-in select_first, or else the user inputting a mask file
	# will result in WDL looking for the literal gs:// URI rather than while the file
	# is localized. Different WDL executors localize files to different places, so the
	# following workaround, while goofy, seems to be the most robust.
	if [[ "~{tbmf}" = "" ]]
	then
		mask="/mask/R00000039_repregions.bed"
	else
		mask="~{tbmf}"
	fi
	
	echo "Pulling diff script..."
	wget https://raw.githubusercontent.com/aofarrel/parsevcf/1.3.1/vcf_to_diff_script.py
	echo "Running script..."
	python3 vcf_to_diff_script.py -v ~{vcf} \
	-d . \
	-tbmf ${mask} \
	-bed ~{bedgraph} \
	-cd ~{min_coverage_per_site}
	
	# if the sample has too many low coverage sites, throw it out
	this_files_info=$(awk -v file_to_check="~{basename_vcf}.diff" '$1 == file_to_check' "~{basename_vcf}.report")
	if [[ ! "$this_files_info" = "" ]]
	then
		# okay, we have information about this file. is it above the removal threshold?
		echo "$this_files_info" > temp
		amount_low_coverage=$(cut -f2 temp)
		
		# account for very tiny numbers
		if [[ "$amount_low_coverage" == *"e"* ]]
		then
			echo "Scientific notation detected, so it's likely this sample is very much passing."
			echo "PASS" >> ERROR
			exit 0
		fi
		
		percent_low_coverage=$(echo "$amount_low_coverage"*100 | bc)
		echo "$percent_low_coverage percent of ~{basename_vcf} is below ~{min_coverage_per_site}x coverage."
		
		# piping an inequality to `bc` will return 0 if false, 1 if true
		is_bigger=$(echo "$amount_low_coverage>~{max_ratio_low_coverage_sites_per_sample}" | bc)
		if [[ $is_bigger == 0 ]]
		then
			# amount of low coverage is BELOW the removal threshold: sample passes
			echo "PASS" >> ERROR
	
		else
			# amount of low coverage is ABOVE the removal threshold: sample fails
			if [[ "~{force_diff}" == "false" ]]
			then
				rm "~{basename_vcf}.diff"
			fi
			pretty_percent=$(printf "%0.2f" "$percent_low_coverage")
			echo FAILURE - "$pretty_percent""%" is above "~{max_ratio_low_coverage_sites_per_sample}""%" cutoff
			echo VCF2DIFF_"$pretty_percent"_PCT_BELOW_"~{min_coverage_per_site}"x_COVERAGE >> ERROR
		fi
	fi
	>>>
	
	runtime {
		cpu: cpu
		docker: "ashedpotatoes/sranwrp:1.1.15"
		disks: "local-disk " + finalDiskSize + " HDD"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	meta {
		author: "Lily Karim (WDLization by Ash O'Farrell)"
	}

	output {
		File? diff = basename_vcf+".diff"
		File? report = basename_vcf+".report"
		String errorcode = read_string("ERROR")
	}
}

task make_diff_legacy {
	input {
		File vcf
		File tbmf
		File cf
		Int cd = 10

		# runtime attributes
		Int addldisk = 10
		Int cpu	= 8
		Int retries	= 1
		Int memory = 16
		Int preempt	= 1
	}
	# estimate disk size
	String basename = basename(vcf)
	Int finalDiskSize = 2*ceil(size(vcf, "GB")) + addldisk

	command <<<
		set -eux pipefail
		mkdir outs
		wget https://raw.githubusercontent.com/lilymaryam/parsevcf/1.0.4/vcf_to_diff_script.py
		python3.10 vcf_to_diff_script.py -v ~{vcf} -d ./outs/ -tbmf ~{tbmf} -cf ~{cf} -cd ~{cd}
		ls -lha outs/
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + finalDiskSize + " SSD"
		docker: "ashedpotatoes/sranwrp:1.1.6"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	meta {
		author: "Lily Karim (WDLization by Ash O'Farrell)"
	}

	output {
		File diff = "outs/"+basename+".diff"
		File report = "outs/"+basename+".report"
	}
}
