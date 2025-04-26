# hox_gene_family
Detect the hox cluster on the grasshopper genome Locusta migratoria

# 1. Detect the Hox gene cluster
## 1.1 Using the genome annotation
   After re-annotated the genome, we found the seven candiates hox family memebers. It's very hard to identfy the hox genes, causing the hox shared a very conserved region but the un-conserved region are variable. So it can't clarify the hox members base on the ortholog method.

## 1.2 tblastx
   Firstly, we take the protein and annotation of hox family memebers from the fruit fly genome. We utilizing the highly conserved region search the corresponding hox candiates in Locusta migratoria. And then the Schistocerca species genome and annotation were used to check the hox candates again by tblastx.
   
   The results are in https://docs.google.com/spreadsheets/d/1kY4OBmMStM12mcQGZ1Ir3BdIPqlKyHbl8GAv9-C__vQ/edit?usp=sharing

### 1.2.1 Download the sequences
wget the Schistocerca species genome from PRJNA772266 and fruit fly Drosophila melanogaster and grep the hox genes name from gff files

          for i in $(cat hox_fullname.genelist) 
          do
          grep "$i" GCF_023864345.2_iqSchSeri2.2_genomic.gff|awk -vFS='\t' '$3=="CDS"{print$1,$9}' >> hox_gene_fullname_iqSchSeri2.2.chr.list
          grep "$i" GCF_023897955.1_iqSchGreg1.2_genomic.gff|awk -vFS='\t' '$3=="CDS"{print$1,$9}' >> hox_gene_fullname_iqSchGreg1.2.chr.list
          grep "$i" GCF_023898315.1_iqSchNite1.1_genomic.gff|awk -vFS='\t' '$3=="CDS"{print$1,$9}' >> hox_gene_fullname_iqSchNite1.1.chr.list
          done

          cat hox_fullname.genelist
            labial
            proboscipedia
            zerknüllt
            deformed
            sex\ combs\ reduced
            fushi tarazu
            antennapedia
            ultrabithorax
            abdominal-A
            abdominal-B
            abdominal
            zerkn
            
          cat get_hox.sh
              for i in $(cat hox.genelist)
              do
              less GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz|grep -w "$i"|awk -vFS='\t' '{print$i,$9}'|grep 'ID=gene' >> hox_gene_fruitfly.chr.list
              done
              
  mutiplt alignment search the hox conserved region
           mafft --auto --clustalout hox_genename_all_chr11.list.pep.fa > hox_genename_all_chr11.list.pep.mafft
          
### 1.2.2tblastx
Take out the core region protein sequences and tblastx on the Locusta migratoria chr11
         makeblastdb -in chr11.fasta  -dbtype nucl -out chr11_nul
         makeblastdb -in chr11.fasta  -dbtype prot -out chr11_pro
         blastn -db chr11_nul -query hox_conserve_region.pep -outfmt 6 -out hox_conserve_region2chr11_pro_outform6
filter the results
      awk -vOFS="\t" '$3>90 {print$2,$9,$10,$1"_"$3}' hox_conserve_region2chr11_pro_outform6 > hox_conserve_region2chr11_pro_outform6.bed
      awk '{if ($2 > $3) print $1"\t"$3"\t"$2"\t"$4; else print $1"\t"$2"\t"$3"\t"$4}' hox_conserve_region2chr11_pro_outform6.bed > hox_conserve_region2chr11_pro_outform6.bed.1
      sed -i 's/CM048754.1/chr11/g' hox_conserve_region2chr11_pro_outform6.bed.1
      bedtools intersect -a hox_conserve_region2chr11_pro_outform6.bed.1 -b ../../locust_chr11.EVM.gtf -wa -wb > hox_conserve_region2chr11_pro_outform6.bed.gtf
      awk '$7=="gene"{print$0}' hox_conserve_region2chr11_pro_outform6.bed.gtf|les

      for i in $(cat hox_ident90_gene.list)
      do
      grep -w $i hox_conserve_region2chr11_pro_outform6.bed.gtf >> $i.list
      done
Finally got the location and gene name of hox genes on the Locusta migratoria chr11.

## 1.3 plot the syteny
### 1.3.1 Plotsr  https://github.com/schneebergerlab/plotsr?tab=readme-ov-file
Using the McscanX results and filtered the 100% and <90% identity alignments.

      #in /home/shangao/Scratch/grasshopper/genome/chr11_maker/chr11tochr11_mcscanx/mcscanx/test
      awk -vOFS="\t" '$3>90 && $3!=100{print$0}' ../xyz.blast > xyz.blast
      
      /home/shangao/software/MCScanX/MCScanX xyz 
      /home/shangao/software/NGenomeSyn/bin/MCScanX2Link.pl ../A.gff ../B.gff  xyz.collinearity 123
      
      awk -vOFS="\t" '{print$0,"SYN"}' 123.A2B.link > ~/Scratch/grasshopper/genome/chr11_maker/hox_family/genomes/plot_genome/chr11_self90.bp
      
      #./
      sed -i 's/chr2/CM048754.1/g' chr11_self90.bp
      sed -i 's/chr1/CM048754.1/g' chr11_self90.bp
      #awk -vOFS="\t" '{if ((2000000 < $2 && $2 < 12000000) || (70000000<$3 && $3<90000000)) print $0,"cl:pink;lw:2;z:2";else print$0,"cl:blue;lw:2;z:1"}' chr11_self93.bp > chr11_self93.f.bp
      python ~/script/python/fliter_plotsr_align.py -s chr11_self90.bp -o chr11_self90.f.bp
      python /home/shangao/software/plotsr/bin/plotsr --bp chr11_self90.f.bp --genomes genomes3.txt -o output_plot1.pdf --markers hox_all.sort.order.bed.4 -W 10 -f 10 -H 3

### 1.3.2 SVbyEye https://htmlpreview.github.io/?https://github.com/daewoooo/SVbyEye/blob/master/man/doc/SVbyEye.html#generate-all-versus-all-minimap-alignments
Easy to use, from minimap2 to plot.
      
      /home/shangao/software/cactus-bin-v2.6.7/bin/minimap2 -x asm20 -c -eqx -secondary=no fruitfly_NT_033777.3.fasta iqSchSeri2.2_NC_064649.1.fasta > fruitfly2iqSchSeri2.2.align
      R
      library(SVbyEye)
      paf.table <- readPaf(paf.file = paf.file,include.paf.tags = TRUE, restrict.paf.tags = "cg")
      filterPaf(paf.table = paf.table, min.align.len = 100000)
      plotMiro(paf.table = paf.table, color.by = "identity")
      dev.off()
# 2. Different expression genes in eggs
## 2.1 Clean and mapping
      Clean
      trim_galore -j 40 -q 30 --fastqc --paired --output_dir ./ /home/zhangtingting/grasshopper/02RNA-Seq/20220721/zen/dsGFP/${i}_1.fq.gz /home/zhangtingting/grasshopper/02RNA-Seq/20220721/zen/dsGFP/${i}_2.fq.gz
      Mapping
      /home/zhangtingting/software/STAR-2.7.11b/source/STAR --genomeDir /home/zhangtingting/grasshopper/01genome/LG_STAR --runThreadN 20 --readFilesIn /home/zhangtingting/grasshopper/02RNA-Seq/20220721/zen/dsLmzen/${i}_1.fq.gz,/home/zhangtingting/grasshopper/02RNA-Seq/20220721/zen/dsLmzen/${i}_2.fq.gz --readFilesCommand zcat --outFileNamePrefix zen_$i --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFilterMismatchNmax 3 --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNoverReadLmax 0.5
## 2.2 Quantity
      for i in 2020_I2N5-dsGFP2_L3_135135Aligned.sortedByCoord.out.bam \
      2020_I2N5-dsGFP3_L3_136136Aligned.sortedByCoord.out.bam \
      2020_I2N5-dszen1_L3_137137Aligned.sortedByCoord.out.bam \
      2020_I2N5-dszen3_L3_138138Aligned.sortedByCoord.out.bam \
      2020_W5d2IN-1_L1_312312Aligned.sortedByCoord.out.bam \
      2d-zen_GFP1-24hAligned.sortedByCoord.out.bam \
      2d-zen_GFP2-24hAligned.sortedByCoord.out.bam \
      2d-zen_GFP3-24hAligned.sortedByCoord.out.bam \
      2d-zen_zen1-24hAligned.sortedByCoord.out.bam \
      2d-zen_zen2-24hAligned.sortedByCoord.out.bam \
      2d-zen_zen3-24hAligned.sortedByCoord.out.bam \
      zen_GFP5Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_GFP1Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_GFP2Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_GFP3Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_GFP4Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_zen2-1Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_zen2-2Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_zen4Aligned.sortedByCoord.out.bam \
      zen_V350058468_L02_zen3Aligned.sortedByCoord.out.bam
      do
      stringtie -p 10 -G /home/zhangtingting/grasshopper/02RNA-Seq/evidenceModeler/evidence_zxm/EVM.all.gff -o $i.gtf /home/zhangtingting/grasshopper/02RNA-Seq/zenmap2genome/$i
      echo $i.gtf >> gtf_list.txt
      done
      stringtie --merge -p 10 -G /home/zhangtingting/grasshopper/02RNA-Seq/evidenceModeler/evidence_zxm/EVM.all.gff -o stringtie_merged.gtf gtf_list.txt
      
      for i in 2020_D5d2IN-1_L1_314314Aligned.sortedByCoord.out.bam \
      2020_I2N5-dsGFP2_L3_135135Aligned.sortedByCoord.out.bam \
      2020_I2N5-dsGFP3_L3_136136Aligned.sortedByCoord.out.bam \
      2020_I2N5-dszen1_L3_137137Aligned.sortedByCoord.out.bam \
      2020_I2N5-dszen3_L3_138138Aligned.sortedByCoord.out.bam \
      2020_W5d2IN-1_L1_312312Aligned.sortedByCoord.out.bam \
      2d-zen_GFP1-24hAligned.sortedByCoord.out.bam \
      2d-zen_GFP2-24hAligned.sortedByCoord.out.bam \
      2d-zen_GFP3-24hAligned.sortedByCoord.out.bam \
      2d-zen_zen1-24hAligned.sortedByCoord.out.bam \
      2d-zen_zen2-24hAligned.sortedByCoord.out.bam \
      2d-zen_zen3-24hAligned.sortedByCoord.out.bam \
      zen_GFP5Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_GFP1Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_GFP2Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_GFP3Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_GFP4Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_zen2-1Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_zen2-2Aligned.sortedByCoord.out.bam \
      zen_V350058468_L01_zen4Aligned.sortedByCoord.out.bam \
      zen_V350058468_L02_zen3Aligned.sortedByCoord.out.bam
      do
      stringtie -e -B -p 30 -G stringtie_merged.gtf -o ballgown/$i.gtf /home/zhangtingting/grasshopper/02RNA-Seq/zenmap2genome/$i
      echo $i.gtf >> ballgown/gtf_list.txt
      done
      sed s/\.gtf//g ballgown/gtf_list.txt > ballgown/gtf_list.txt1
      paste -d '\t' ballgown/gtf_list.txt1 ballgown/gtf_list.txt > ballgown/gtf_list.txt2
      cd ballgown
      python prepDE.py3 -i gtf_list.txt2

## 2.3 Plot of DEGs
      library("DESeq2")
      countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
      setwd("D:/R_space/grasshopper_eggs_deg")
      colData <- read.table('sample.list.txt', header=T,row.names=1)
        

      
