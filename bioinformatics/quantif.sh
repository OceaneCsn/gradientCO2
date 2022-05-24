cd /data/ocassan_data/gradientCO2/bams/
for i in `ls *.bam`
do
    echo $i
    name=${i%%Aligned.sortedByCoord.out.bam} ###Vire le .bam
    echo $name
    htseq-count -f bam -r pos $i /data/ocassan_data/TAIR10/TAIR10_GFF3_genes.gff --idattr=ID --type=gene --stranded=no > ../quantif_unstranded_genes/$name.txt &
done 
