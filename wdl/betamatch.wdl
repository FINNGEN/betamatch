version 1.0

task match_betas {

    input {
        Array[File] ref_files
        Array[File] summary_files
        Array[File] tbi_indexes
        String docker

        Array[String] column_names_ref
        Array[String] column_names_other
        Float pval_threshold
        String zones
        String xlabel
        String ylabel

        String out_f = "out_f"
        Int disk_size = ceil((size(ref_files, "GB") + size(summary_files, "GB")) * 1.2) + 5
    }

    command <<<
        
        set -euxo pipefail

        paste ~{write_lines(ref_files)} ~{write_lines(summary_files)} > matchfile
        mkdir ~{out_f}
        
        betamatch.py --info-ext ~{sep=" " column_names_ref} --info-fg ~{sep=" " column_names_other} --match-file matchfile --output-folder ~{out_f} --pval-filter ~{pval_threshold}
        corrplot.py ~{out_f} --fields unif_beta_ext unif_beta_fg --se-fields ~{column_names_ref[6]}_ext ~{column_names_ref[6]}_fg --x-title "~{xlabel}" --y-title "~{ylabel}" --pval_field ~{column_names_ref[5]}_ext --pval_threshold ~{pval_threshold} --out "output.pdf"

    >>>

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "6 GB"
        disks: "local-disk ~{disk_size} HDD"
        zones: "~{zones}"
        preemptible: 2
    }

    output {
        Array[File] out = glob("out_f/*.betas.tsv")
        File corrplot = "output.pdf"
        File r2_table = "r2_table.tsv"
    }
}

workflow betamatch {

    input {
        String docker
        File match_file

        Array[Array[String]] files = read_tsv(match_file)
        Array[File] ref_files = transpose(files)[0]
        Array[File] summary_files = transpose(files)[1]
    }
    
    scatter (s in range(length(summary_files))) {#scatter magic, it just works
        String t = summary_files[s] + ".tbi"
    }

    call match_betas {
        input:
            ref_files = ref_files,
            summary_files = summary_files,
            tbi_indexes = t,
            docker = docker
    }

    output {
        Array[File] betas = match_betas.out
        File corrplot = match_betas.corrplot
        File r2_table = match_betas.r2_table
    }
}