version 1.0

task make_mask_and_diff_and_process_metadata {
	# For basic information, see make_mask_and_diff.
	#
	# This version additionally processes some metadata fields. If the sample is passing,
	# these fields will be the sample value they were passed in. If the sample is failing,
	# these fields will be undefined. This is a workaround for upstream tasks/pipelines
	# that take in more than one sample at a time and its metadata.
	#
	# "Why can't we just write a TSV file?" Because CalTBNet needs to run per-sample due to
	# how the data table is set up and due to Terra limitations, so this task runs per-sample,
	# so you'd have one two-line TSV file per sample, and then you'd have to localize all of 
	# those thousands of TSV files into a downstream task, and then concatenate them.
	#
	# There are alternatives to doing this, so this is currently not in the main version of myco.
	input {
		File bam
		Boolean force_diff = false
		Boolean histograms = false
		Float max_ratio_low_coverage_sites_per_sample # for sra, default was 0.05
		Int min_coverage_per_site
		File? tbmf
		File vcf

		# metadata key-value pairs; key is the name of the field, value is... value
		# there are complex types that could in theory do this, but Cromwell's
		# handling of optional values in complex types is buggy, so it's much safer
		# to quasi-hardcode the number of metadata fields like this, even if it's cringe
		String? a_key
		String? a_value
		String? b_key
		String? b_value
		String? c_key
		String? c_value
		String? d_key
		String? d_value
		String? e_key
		String? e_value
		String? f_key
		String? f_value

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
		
		# account for very tiny numbers (very big numbers should be impossible)
		if [[ "$amount_low_coverage" == *"e"* ]]
		then
			echo "Scientific notation detected, so it's likely this sample is very much passing."
			echo "PASS" >> ERROR
		else
			percent_low_coverage=$(echo "$amount_low_coverage"*100 | bc)
			maximium_percent_low_coverage=$(echo "~{max_ratio_low_coverage_sites_per_sample}*100" | bc)
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
				echo FAILURE - "$pretty_percent""%" is above "~{min_coverage_per_site}""%" cutoff
				echo VCF2DIFF_"$pretty_percent"_PCT_BELOW_~{min_coverage_per_site}x_COVERAGE_\(MAX_"$maximium_percent_low_coverage"_PCT\) >> ERROR

				end=$(date +%s)
				seconds=$(echo "$end - $start" | bc)
				minutes=$(echo "$seconds" / 60 | bc)
				echo "Finished in about $minutes minutes ($seconds sec)) -- although we QC failed, we'll still write metadata"
			fi
		fi
	fi

	# this executes even if a sample is failing
	python3 << CODE
	a_key =  "~{a_key}"
	a_value = "~{a_value}"
	b_key =  "~{b_key}"
	b_value = "~{b_value}"
	c_key =  "~{c_key}"
	c_value = "~{c_value}"
	d_key =  "~{d_key}"
	d_value = "~{d_value}"
	e_key =  "~{e_key}"
	e_value = "~{e_value}"
	f_key = "~{f_key}"
	f_value = "~{f_value}"

	with open('ERROR') as f:
		status = f.readline()

	valid_keys = []
	for key in [a_key, b_key, c_key, d_key, e_key, f_key]:
		if key == '' or key == ' ' or key == "'":
			print(f"key [{key}] effectively is undefined")
			key = "UNDEFINED"
		print(f"adding {key} to valid_keys")
		valid_keys.append(key.strip("'").strip('"'))
	print(f"valid_keys: {valid_keys}")
	valid_values = []
	for value in [a_value, b_value, c_value, d_value, e_value, f_value]:
		if value == '' or value == ' ' or value == "'":
			print(f"value [{value}] effectively is undefined")
			value = "UNDEFINED"
		print(f"adding {value} to valid_values")
		valid_values.append(value.strip("'").strip('"'))
	print(f"valid_values: {valid_values}")

	assert len(valid_keys) == len(valid_values)
	metadata_dict = dict(zip(valid_keys, valid_values))
	valid_metadata_dict = dict()
	for keys, values in metadata_dict.items():
		if keys == "UNDEFINED" and values == "UNDEFINED":
			print("Keys and values is undefined, dropping")
			continue
		elif keys == "UNDEFINED": # and values does not
			print(f"WARNING: Got metadata value {value} with undefined key")
			continue
		else:
			# it's okay if value is undefined
			print(f"{keys}: {values}")
			valid_metadata_dict[keys] = values
	
	# turn this into something WDL can use
	# for maximum compatibility, we're going to try tabs AND commas
	
	header_tsv = "sample\tstatus\t" + "\t".join(valid_metadata_dict.keys())
	body_tsv = f"~{basename_vcf}\t{status}\t" + "\t".join(valid_metadata_dict.values())
	metadata_tsv = header_tsv + "\n" + body_tsv
	header_tsv, body_tsv, metadata_tsv = header_tsv[:-1], body_tsv[:-1], metadata_tsv[:-1]

	header_csv = "sample,status," + ",".join(valid_metadata_dict.keys())
	body_csv = f"~{basename_vcf},{status}," + ",".join(valid_metadata_dict.values())
	metadata_csv = header_csv + "\n" + body_csv
	header_csv, body_csv, metadata_csv = header_csv[:-1], body_csv[:-1], metadata_csv[:-1]

	for string, file in {header_tsv: "header.tsv", body_tsv: "body.tsv", metadata_tsv: "metadata.tsv", header_csv: "header.csv", body_csv: "body.csv", metadata_csv: "metadata.csv"}.items():
		with open(file, 'w') as f:
			f.write(string)

	CODE

	# how long did this take?
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
		String meta_header_tsv = read_string("header.tsv")
		String meta_header_csv = read_string("header.csv")
		String meta_values_tsv = read_string("body.tsv")
		String meta_values_csv = read_string("body.csv")
		String meta_full_tsv = read_string("metadata.tsv")
		String meta_full_csv = read_string("metadata.csv")
		String errorcode = read_string("ERROR")
	}
}

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
		
		# account for very tiny numbers (very big numbers should be impossible)
		if [[ "$amount_low_coverage" == *"e"* ]]
		then
			echo "Scientific notation detected, so it's likely this sample is very much passing."
			echo "PASS" >> ERROR
		else
			percent_low_coverage=$(echo "$amount_low_coverage"*100 | bc)
			maximium_percent_low_coverage=$(echo "~{max_ratio_low_coverage_sites_per_sample}*100" | bc)
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
				echo FAILURE - "$pretty_percent""%" is above "~{min_coverage_per_site}""%"
				echo VCF2DIFF_"$pretty_percent"_PCT_BELOW_~{min_coverage_per_site}x_COVERAGE_\(MAX_"$maximium_percent_low_coverage"_PCT\) >> ERROR
			fi
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
			echo VCF2DIFF_"$pretty_percent"_PCT_BELOW_~{min_coverage_per_site}x_COVERAGE >> ERROR
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
