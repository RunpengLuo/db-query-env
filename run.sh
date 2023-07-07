#!/bin/bash


base="/Users/luorunpeng/Downloads/all-e/Research/project-database_query/db_query"

ref1=${base}/data/comb1_YEE_vector.fasta
ref2=${base}/data/bbr_comb23_promoter.fasta
ref3=${base}/data/bbr_comb46_signal_peptide.fasta
ref4=${base}/data/bbr_comb1_hsa.fasta 
ref5=${base}/data/comb1_xaMF.fasta
ref6=${base}/data/bbr_comb6_terminator.fasta

qry=${base}/20230504.Plasmid.Mixture1.BC81.fastq

# mkdir k21_test_new
# cd k21_test_new
python ${base}/main.py build -fa 21 k21_db_new $ref1 $ref2 $ref3 $ref4 $ref5 $ref6
python ${base}/main.py query -fq k21_db_new k21_query_csv $qry
# python ${base}/main.py simplify "comps_idx.txt" "query_res.txt"


# k-mer size: 21,33,55,77,99,127
# mkdir k21_test
# cd k21_test
# python ${base}/main.py build -fa 21 $ref1 $ref2 $ref3 $ref4 $ref5 $ref6
# python ${base}/main.py query -fq "comps_db.txt" $qry
# python ${base}/main.py simplify "comps_idx.txt" "query_res.txt"

# cd ..

# mkdir k33_test
# cd k33_test
# python ${base}/main.py build -fa 33 $ref1 $ref2 $ref3 $ref4 $ref5 $ref6
# python ${base}/main.py query -fq "comps_db.txt" $qry
# python ${base}/main.py simplify "comps_idx.txt" "query_res.txt"

# cd ..
# mkdir k36_test
# cd k36_test
# python ${base}/main.py build -fa 36 $ref1 $ref2 $ref3 $ref4 $ref5 $ref6
# python ${base}/main.py query -fq "comps_db.txt" $qry
# python ${base}/main.py simplify "comps_idx.txt" "query_res.txt"

# cd ..