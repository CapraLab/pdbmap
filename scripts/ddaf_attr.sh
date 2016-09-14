# Generates a deltaDAF attribute file from two DAF attribute files
#$1 Structure ID
#$2 Biounit
#$3 Reference Population
#$4 Alternate Population
paste $1_biounit$2_vars_$3_daf.attr $1_biounit$2_vars_$4_daf.attr | awk -v OFS="\t" '{if (NR==1) print $1" "$2" "$3" "$4"_"$8; else if (NR==2) print $1" "$2"_"$4; else if (NR==3) print $1" "$2" "$3; else if (NR<5) print $1" "$2; else print "",$1,$2-$4}' > $1_biounit$2_vars_$3_$4_ddaf.attr
