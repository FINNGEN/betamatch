task match_betas{
    Array[Array[String]] match_file
    Array[File] tbi_indexes
    Array[String] ext_files = transpose(match_file)[0]
    Array[File] summary_stat_files = transpose(match_file)[1] 
    Array[String] column_names_ext
    Array[String] column_names_fg
    String docker
    String out_f = "out_f"
    String ext_repo_url
    String ext_repo_branch
    Float pval_threshold
    command <<<
        #download github repo to ext_repo
        git clone -b ${ext_repo_branch} ${ext_repo_url} ext_repo
        #add ext_repo/data/ to the beginning of ext_files
        cat ${write_lines(ext_files)}|sed -e "s/^/ext_repo\/data\//g" > exts
        #combine the different files into one match file
        paste exts ${write_lines(summary_stat_files)} > matchfile
        mkdir ${out_f}
        
        betamatch.py --info-ext ${sep=" " column_names_ext} --info-fg ${sep=" " column_names_fg} --match-file matchfile --output-folder ${out_f} --pval-filter ${pval_threshold}
        corrplot.py ${out_f} --fields unif_beta_fg unif_beta_ext --se-fields se_fg se_ext --x-title "FinnGen beta" --y-title "External beta" --pval_field pval_ext --pval_threshold ${pval_threshold} --out "output.pdf"
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
    String ext_repo_url
    String ext_repo_branch
    Array[File] temp = transpose(files)[1] 
    Float pval_threshold
    Array[String] column_names_ext
    Array[String] column_names_fg
    scatter (s in range( length( temp) ) ){#scatter magic, it just works
        String t = sub(temp[s],".gz",".gz.tbi") 
    }
    call match_betas {
        input: match_file = files, docker=docker, tbi_indexes=t, pval_threshold = pval_threshold, ext_repo_url = ext_repo_url, ext_repo_branch = ext_repo_branch, column_names_ext = column_names_ext, column_names_fg = column_names_fg
    }
}