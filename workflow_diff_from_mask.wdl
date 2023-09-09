version 1.0

import "./vcf_to_diff.wdl" as tasks

workflow Diff_From_VCF_And_Mask {
    input {
        # index 1 of vcf must correlate with index 1 of bedgraph, etc
        Array[File] vcfs
        Array[File] bedgraphs
        File? tbmf
    }
    
    scatter(vcfs_and_beds in zip(vcfs, bedgraphs)) {
        call tasks.make_diff_from_vcf_and_mask as diffmaker9000 {
            input:
                vcf = vcfs_and_beds.left,
                bedgraph = vcfs_and_beds.right,
                tbmf = tbmf,
                force_diff = true
        }
    }
    
    output {
        Array[File?] diffs = diffmaker9000.diff
        Array[File?] reports = diffmaker9000.report
    }
}

