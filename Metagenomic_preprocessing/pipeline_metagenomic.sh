#### Pipeline for metagenomic data processing
Human_DB=./database/GRCh38/G38 # human reference database
silva_db_kraken2=./database/silva_db_kraken2
ncbi_db_kraken2=./database/ncbi_db_kraken2
repair_sh=./tools/repair.sh

## interface
ori_data_dir=$1
ori_data_suffix1=$2
ori_data_suffix2=$3
metalistfile=$4

## some dirs to save temp files
tmp_out=$ori_data_dir/tmpfiles
clean_out=$ori_data_dir/clean_reads
raw_dir=$ori_data_dir/raw_seqs
mkdir $tmp_out $clean_out $raw_dir

## some parameters for kraken2
silva_db_confidence=0.05
ncbi_db_confidence=0 # GUT -- 0.05 : refer to: Benchmarking Metagenomics Tools for Taxonomic Classification

t_ge_num_ncbi=10 # gut -- 100; soil -- 10
t_ge_num_silva=10

## start pipeline
for demo_nm in `cat $metalistfile`
do
    out_dir_this_srr=$clean_out'/'$demo_nm
    mkdir -p $out_dir_this_srr

    fq1_ori=$ori_data_dir'/'$demo_nm$ori_data_suffix1
    fq2_ori=$ori_data_dir'/'$demo_nm$ori_data_suffix2

    fq1=$raw_dir/$demo_nm'_1.fastq'
    fq2=$raw_dir/$demo_nm'_2.fastq'

    if [ ! -f $out_dir_this_srr/$demo_nm'_step4_ncbi_rp.txt' ]; then # $tmp_out/$demo_nm'_step2_hostRemoved.1.fastq'
        echo "Processed file not found: start from begining!"
        rm $tmp_out/$demo_nm'_step1'* # rm tmpfiles if exists

        if [ "${fq1_ori##*.}" == "gz" ]; then
            echo "${fq1_ori##*.}"
            # Decompose the gz files.
            echo "Decompose the gz files!"
            gunzip -c $fq1_ori > $fq1
            gunzip -c $fq2_ori > $fq2
        else
            echo "Parse the sra file and move it to raw dir!"
            fasterq-dump --split-files $ori_data_dir'/'$demo_nm'/' -O $ori_data_dir'/'

            mv $fq1_ori $fq1
            mv $fq2_ori $fq2
        fi

        # Normalize the reads header to /1 and /2 using awk.
        awk '{if(NR % 4 == 1 || NR % 4 == 3){if($0 !~ /\/1$/){sub("2:N:0:", "1:N:0:"); split($0, fields, " "); $0 = fields[1] "/1"; print}else{print}}else{print}}' $fq1 > ${fq1}.temp
        awk '{if(NR % 4 == 1 || NR % 4 == 3){if($0 !~ /\/2$/){sub("2:N:0:", "1:N:0:"); split($0, fields, " "); $0 = fields[1] "/2"; print}else{print}}else{print}}' $fq2 > ${fq2}.temp
        mv "${fq1}.temp" "$fq1"
        mv "${fq2}.temp" "$fq2"
        # Use repair.sh to solve the unmathed reads in -1 and -2 files.
        $repair_sh in1=$fq1 in2=$fq2 out=${fq1}.temp out2=${fq2}.temp outs=${fq1}.single.fq
        mv "${fq1}.temp" "$fq1"
        mv "${fq2}.temp" "$fq2"

        ## 1. QC, primer-header
        # Use fastp to do QC and remove primer-header, the deatils of parameters can be found by 'fastp -h' or in github readme.
        # Note: the adapter seq will be automatically detected, but you can also specify it by yourself.
        # length: 100
        fastp -i $raw_dir/$demo_nm'_1.fastq' -I $raw_dir/$demo_nm'_2.fastq' -o $tmp_out/$demo_nm'_step1_1.fastq' -O $tmp_out/$demo_nm'_step1_2.fastq' \
            -h $tmp_out/$demo_nm'_step1.html' \
            -j $tmp_out/$demo_nm'_step1.json' \
            --correction --overlap_len_require 20 --overlap_diff_limit 5 --overlap_diff_percent_limit 20 \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 30 \
            --length_required 85 \
            --length_limit 500 \
            --complexity_threshold 30 \
            --detect_adapter_for_pe \
            --cut_tail --cut_tail_mean_quality 15 \
            --dont_overwrite \
            --low_complexity_filter \
            -w 8
            # --merge \
            # --merged_out $tmp_out/$demo_nm'_step1_merged.fastq' \
        rm $raw_dir/$demo_nm'_1.fastq' $raw_dir/$demo_nm'_2.fastq'

        ## 2. Remove host seqs
        # Use bowtie2 to remove host seqs, the deatils of parameters can be found by 'bowtie2 -h' or in github readme.
        bowtie2 -x $Human_DB -1 $tmp_out/$demo_nm'_step1_1.fastq' -2 $tmp_out/$demo_nm'_step1_2.fastq' \
            --un-conc $tmp_out/$demo_nm'_step2_hostRemoved.fastq' \
            --very-fast \
            -p 16 \
            --very-fast-local --quiet > $tmp_out/$demo_nm'_step2_ot.txt'

        ## 3. Kraken2 to classify
        kraken2 --paired --db $silva_db_kraken2 \
            --unclassified-out $out_dir_this_srr/$demo_nm'_step3_unclassified#.fq'  \
            --classified-out $out_dir_this_srr/$demo_nm'_step3_classified#.fq' \
            --report $out_dir_this_srr/$demo_nm'_step3_rp.txt' --output $out_dir_this_srr/$demo_nm'_step3_ot.txt' \
            --threads 16 \
            --confidence $silva_db_confidence \
            $tmp_out/$demo_nm'_step2_hostRemoved.1.fastq' $tmp_out/$demo_nm'_step2_hostRemoved.2.fastq'

        kraken2 --paired --db $ncbi_db_kraken2 \
            --unclassified-out $out_dir_this_srr/$demo_nm'_step3_ncbi_unclassified#.fq'  \
            --classified-out $out_dir_this_srr/$demo_nm'_step3_ncbi_classified#.fq' \
            --report $out_dir_this_srr/$demo_nm'_step3_ncbi_rp.txt' --output $out_dir_this_srr/$demo_nm'_step3_ncbi_ot.txt' \
            --threads 16 \
            --confidence $ncbi_db_confidence \
            $tmp_out/$demo_nm'_step2_hostRemoved.1.fastq' $tmp_out/$demo_nm'_step2_hostRemoved.2.fastq'
    
        ## 3.5. extract seqs classified to Bacteria and do it again.
        # extract silva classified seqs
        extract_kraken_reads.py -k $out_dir_this_srr/$demo_nm'_step3_ot.txt' -s1 $out_dir_this_srr/$demo_nm'_step3_classified_1.fq' -s2 $out_dir_this_srr/$demo_nm'_step3_classified_2.fq' \
            -o $out_dir_this_srr/$demo_nm'_extraBacteria_1.fq' -o2 $out_dir_this_srr/$demo_nm'_extraBacteria_2.fq' \
            -r $out_dir_this_srr/$demo_nm'_step3_rp.txt' \
            -t 3 --include-children # extract bacteria reads
        # extract ncbi classified seqs
        extract_kraken_reads.py -k $out_dir_this_srr/$demo_nm'_step3_ncbi_ot.txt' -s1 $out_dir_this_srr/$demo_nm'_step3_ncbi_classified_1.fq' -s2 $out_dir_this_srr/$demo_nm'_step3_ncbi_classified_2.fq' \
            -o $out_dir_this_srr/$demo_nm'_extraBacteria_ncbi_1.fq' -o2 $out_dir_this_srr/$demo_nm'_extraBacteria_ncbi_2.fq' \
            -r $out_dir_this_srr/$demo_nm'_step3_ncbi_rp.txt' \
            -t 2 --include-children # extract bacteria reads, here bacteria taxid is 2
        # classify the extraBacteria seqs
        kraken2 --paired --db $silva_db_kraken2 \
            --unclassified-out $out_dir_this_srr/$demo_nm'_step3.5_ExtraBac_unclassified#.fq'  \
            --classified-out $out_dir_this_srr/$demo_nm'_step3.5_ExtraBac_classified#.fq' \
            --report $out_dir_this_srr/$demo_nm'_step3.5_ExtraBac_rp.txt' --output $out_dir_this_srr/$demo_nm'_step3.5_ExtraBac_ot.txt' \
            --threads 16 \
            --confidence $silva_db_confidence \
            $out_dir_this_srr/$demo_nm'_extraBacteria_1.fq' $out_dir_this_srr/$demo_nm'_extraBacteria_2.fq' # only classify the bacteria seqs
        kraken2 --paired --db $ncbi_db_kraken2 \
            --unclassified-out $out_dir_this_srr/$demo_nm'_step3.5_ExtraBac_ncbi_unclassified#.fq'  \
            --classified-out $out_dir_this_srr/$demo_nm'_step3.5_ExtraBac_ncbi_classified#.fq' \
            --report $out_dir_this_srr/$demo_nm'_step3.5_ncbi_ExtraBac_rp.txt' --output $out_dir_this_srr/$demo_nm'_step3.5_ncbi_ExtraBac_ot.txt' \
            --threads 16 \
            --confidence $ncbi_db_confidence \
            $out_dir_this_srr/$demo_nm'_extraBacteria_ncbi_1.fq' $out_dir_this_srr/$demo_nm'_extraBacteria_ncbi_2.fq' # only classify the bacteria seqs
    fi
    
    ## 4. Get the relative abun matrix
    bracken -d $silva_db_kraken2 -i $out_dir_this_srr/$demo_nm'_step3.5_ExtraBac_rp.txt' \
        -o $out_dir_this_srr/$demo_nm'_step4_abun.txt' -w $out_dir_this_srr/$demo_nm'_step4_rp.txt' -r 150 \
        -l G -t $t_ge_num_silva # generate relative abun matrix to genus level
    bracken -d $ncbi_db_kraken2 -i $out_dir_this_srr/$demo_nm'_step3.5_ncbi_ExtraBac_rp.txt' \
        -o $out_dir_this_srr/$demo_nm'_step4_ncbi_abun.txt' -w $out_dir_this_srr/$demo_nm'_step4_ncbi_rp.txt' -r 150 \
        -l G -t $t_ge_num_ncbi # generate relative abun matrix to genus level

    ## rm tmp files
    rm $tmp_out/$demo_nm'_step2_ot.txt' $tmp_out/$demo_nm'_step1'*'.fastq' # $tmp_out/$demo_nm*'.fastq'
    rm $tmp_out/$demo_nm'_step2_hostRemoved.1.fastq' $tmp_out/$demo_nm'_step2_hostRemoved.2.fastq'
    rm $out_dir_this_srr/$demo_nm'_step3'*'classified'*'.fq' $out_dir_this_srr/$demo_nm'_extraBacteria'*'.fq' $out_dir_this_srr/$demo_nm'_step3.5_ExtraBac_unclassified'*'.fq'
    rm $out_dir_this_srr/$demo_nm*'_ot.txt'
done