sed 's/\,/\t/g' res_treatment_ins_vs_con_all_genes.csv > outfile

#cut -f7,3,6 outfile > res_treatment_ins_vs_con.bsf

awk '{ print $7 "\t" $3 "\t" $6}' outfile > res_treatment_ins_vs_con.bsf

#open vim remove header

awk  '$2!="NA"' # remove NAs
