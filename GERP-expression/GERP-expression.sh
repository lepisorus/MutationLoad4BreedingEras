#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=96:00:00
#SBATCH --job-name=exp
#SBATCH --output=R_%j.out 
#SBATCH --error=R_%j.err
#SBATCH --mail-user=sunshichao@caas.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module load python
module load perl
module load bcftools
module load vcftools


##Allele frequency calculation
vcftools --vcf CN1.vcf --freq --out CN1.freq
vcftools --vcf CN2.vcf --freq --out CN2.freq
vcftools --vcf CN3.vcf --freq --out CN3.freq

##Obtain MAF
cat <(cut -f1,2,5,6 CN1.freq.frq |sed 's/\:/\t/g' |awk '{if($4<$6) print$0,$7=$4}') <(cut -f1,2,5,6 CN1.freq.frq |sed 's/\:/\t/g' |awk '{if($4>$6) print$0,$7=$6}') |awk '{print$1"_"$2"\t"$7}' >CN1.MAF.txt
cat <(cut -f1,2,5,6 CN2.freq.frq |sed 's/\:/\t/g' |awk '{if($4<$6) print$0,$7=$4}') <(cut -f1,2,5,6 CN2.freq.frq |sed 's/\:/\t/g' |awk '{if($4>$6) print$0,$7=$6}') |awk '{print$1"_"$2"\t"$7}' >CN2.MAF.txt
cat <(cut -f1,2,5,6 CN3.freq.frq |sed 's/\:/\t/g' |awk '{if($4<$6) print$0,$7=$4}') <(cut -f1,2,5,6 CN3.freq.frq |sed 's/\:/\t/g' |awk '{if($4>$6) print$0,$7=$6}') |awk '{print$1"_"$2"\t"$7}' >CN3.MAF.txt

##Take the SNP intersection of three populations	
awk '{if($2!=0) print$0}' CN1.MAF.txt|cut -f1 >CN1.MAF.list
awk '{if($2!=0) print$0}' CN2.MAF.txt|cut -f1 >CN2.MAF.list
awk '{if($2!=0) print$0}' CN3.MAF.txt|cut -f1 >CN3.MAF.list

sort CN1.MAF.list CN2.MAF.list |uniq -d >CN1_2.txt
sort CN1_2.txt CN3.MAF.list |uniq -d >CN1_2_3.txt

##Obtain new file of allele frequency
perl ../get.pl CN1.MAF.freq CN1_2_3.txt >CN1.MAF.common.snps.freq
perl ../get.pl CN2.MAF.freq CN1_2_3.txt >CN2.MAF.common.snps.freq
perl ../get.pl CN3.MAF.freq CN1_2_3.txt >CN3.MAF.common.snps.freq

##The bed file of the expression genes TSS-2kb
cat geneexpression.txt |cut -f1 >gene.list
awk '{print$4"\t"$1"\t"$2"\t"$3}' zea_mays.protein_coding.bed >zea_mays.protein_coding.query
perl get.pl zea_mays.protein_coding.query gene.list >gene.expression.sites
awk '{print$2"\t"$3-2000"\t"$3"\t"$1}' gene.expression.sites |sort -k1n,1 -k2n,2 >gene.expression.snps.bed
bedtools intersect -a gene.expression.snps.bed -b all.snps.bed -wa -wb |awk '{print$5"_"$6"\t"$4}' >gene.expression.id.txt
cut -f1 gene.expression.id.txt >SNPs.list

##SNP sites contained at TSS-2Kb of the expression gene
perl get.pl CN1.MAF.txt CN1.MAF.list >CN1.MAF
perl get.pl CN2.MAF.txt CN2.MAF.list >CN2.MAF
perl get.pl CN3.MAF.txt CN3.MAF.list >CN3.MAF
paste CN1.MAF CN2.MAF CN3.MAF |awk '{print$1"\t"$2"\t"$4"\t"$6}' >SNP.MAFCN1_CN2_CN3.txt

#gerp value
le allSNPs.GERP0.list |awk '{print$1"_"$2"\t"$3}' >ID.gerp.txt
le gene.id.name.txt |cut -f1 >query
perl ../get.pl ID.gerp.txt query >gene.exp.gerp.txt
paste gene.id.name.txt gene.exp.gerp.txt | cut -f1,2,4 >ID.GENE.GERP.txt
le ID.GENE.GERP.txt |cut -f2 |sort|uniq >gene.uniq

##The hightest value of GERP 
for i in `cat gene.uniq`;do printf "cat geneid.snp.gerp.txt |grep "$i" |sort -nrk3,3 |head -1 >./tmp/$i.txt\n";done >run.sh
cat *.txt >ID.GENE.GERP.highest.txt
cut -f1 ID.GENE.GERP.highest.txt >ID.GERP.highest.query

perl ../get.pl CN1.MAF ID.GERP.highest.query >CN1.GERP.highest.freq
perl ../get.pl CN2.MAF ID.GERP.highest.query >CN2.GERP.highest.freq
perl ../get.pl CN3.MAF ID.GERP.highest.query >CN3.GERP.highest.freq

##expression
le CN1.expression.txt |grep -v "model" |cut -f2 |awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' >CN1.average.expression.txt
le CN2.expression.txt |grep -v "model" |cut -f2 |awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' >CN2.average.expression.txt
le CN3.expression.txt |grep -v "model" |cut -f2 |awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' >CN3.average.expression.txt

grep -v "model" CN1.expression.txt|cut -f1 >CN1.1.txt
grep -v "model" CN2.expression.txt|cut -f1 >CN1.2.txt
grep -v "model" CN3.expression.txt|cut -f1 >CN1.3.txt
paste CN1.1.txt CN1.average.expression.txt CN2.average.expression.txt CN3.average.expression.txt |head -18521 >geneid.expCN1_CN2_CN3.txt

##Combine geneID and SNPs
cut -f1 SNP.frqCN1_CN2_CN3.txt >snp.list
perl ../get.pl gene.id.name.txt snp.list >query.result
cut -f2 query.result>geneid.query
perl ../get.pl geneid.expCN1_CN2_CN3.txt geneid.query > same.sort.snp_and_geneid.txt
paste SNP.frqCN1_CN2_CN3.txt same.sort.snp_and_geneid.txt >final.plot.r2.txt
less final.plot.r2.txt |awk '{if($6>=1 || $7>=1 || $8>=1) print$0}' |awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' >final.plot.maf.filter.txt