task match_betas{
    Array[Array[String]] match_file
    Array[File] tbi_indexes
    Array[File] ext_files = transpose(match_file)[0]
    Array[File] summary_stat_files = transpose(match_file)[1]
    Array[String] column_names=["\"#chrom\"","pos","ref","alt","beta","pval"]
    String docker
    Int mem
    String out_f = "out_f"
    Float pval_threshold

    command <<<
        #combine the different files into one match file
        paste exts ${write_lines(summary_stat_files)} > matchfile
        mkdir ${out_f}
        betamatch.py --info ${sep=" " column_names} --match-file matchfile --output-folder ${out_f}
        corrplot.py ${out_f} --fields unif_beta_fg unif_beta_ext --se-fields sebeta se --x-title "FinnGen beta" --y-title "External beta" --pval_field pval_ext --pval_threshold ${pval_threshold} --out "output.pdf"
    >>>

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "${mem} GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b"
        preemptible: 2
    }

    output {
        Array[File] out = glob("out_f/*.betas.tsv")
        File corrplot = "output.pdf"
        File r2_table = "r2_table.tsv"
    }
}

workflow betamatch{
    String docker
    Int mem
    File match_file
    Array[Array[String]] files = read_tsv(match_file)
    Array[File] temp = transpose(files)[1]
    Float pval_threshold
    scatter (s in range( length( temp) ) ){#scatter magic, it just works
        String t = sub(temp[s],".gz",".gz.tbi")
    }
    call match_betas {
        input: match_file = files, docker=docker, mem = mem, tbi_indexes=t, pval_threshold = pval_threshold
    }
}
