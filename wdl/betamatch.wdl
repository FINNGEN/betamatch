task match_betas{
    Array[Array[File]] match_file
    Array[File] temp = transpose(match_file)[0] 
    Array[File] tbi_indexes
    Array[String] column_names=["\"#chrom\"","pos","ref","alt","beta","pval"]
    String docker
    String out_f = "out_f"
    Float pval_threshold
    command <<<
        mkdir ${out_f}
        betamatch.py --info ${sep=" " column_names} --match-file ${write_tsv(match_file)} --output-folder ${out_f}
        corrplot.py ${out_f} --fields unif_beta_fg unif_beta_ext --se-fields sebeta se --x-title "FinnGen beta" --y-title "External beta" --pval_field pval_ext --pval_threshold ${pval_threshold} --out "output.pdf"
    >>>

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "3 GB"
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
    File match_file
    Array[Array[String]] files = read_tsv(match_file)
    Array[File] temp = transpose(files)[1] 
    Float pval_threshold
    scatter (s in range( length( temp) ) ){#scatter magic, it just works
        String t = sub(temp[s],".gz",".gz.tbi") 
    }
    call match_betas {
        input: match_file = files, docker=docker, tbi_indexes=t, pval_threshold = pval_threshold
    }
}