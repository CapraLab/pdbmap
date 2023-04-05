# Every time we get new ensembl database, the very huge files must be split to have any sense of progress.
echo WARNING DO NOT CHECK IN WITH MYSQLIMPORT PASSWORD IN COMMAND LINE
version="108_38"
echo ENSEMBL version is $version
echo  # cd to homo_sapiens_variation_108_38 from
cd /dors/capra_lab/data/ensembl/2022-11-17/homo_sapiens_variation_$version
mkdir -p transcript_variation_split
split --verbose -l 5000000 -a 4 transcript_variation.txt transcript_variation_split/transcript_variation.piece.
# Now that we've split, carefully do mysqlimport of each of the split-out files.  Take care with passwords
cd transcript_variation_split
piece_filenames=( `echo transcript_variation.piece.*` )
for transcript_variation_piece in "${piece_filenames[@]}"
do
echo $transcript_variation_piece
mv --verbose $transcript_variation_piece transcript_variation.txt
echo $transcript_variation_piece ": begin mysqlimport"
mysqlimport --fields-terminated-by='\t' --fields-escaped-by=\\ homo_sapiens_variation_$version -preplace_with_pw -L transcript_variation.txt | tee -a mysqlimport.log
if [$? -eq 0 ]; then
    echo "mysqlimport successfully returned 0.  Removing input piece file"
    rm --verbose transcript_variation.txt
else
    echo "mysqlimport failed.  Retaining input piece file"
    mv --verbose transcript_variation.txt $transcript_variation_piece
fi
done
