# hox_gene_family
Detect the hox cluster on the grasshopper genome Locusta migratoria

1. Detect the Hox gene cluster
1.1 Using the genome annotation
   After re-annotated the genome, we found the seven candiates hox family memebers. It's very hard to identfy the hox genes, causing the hox shared a very conserved region but the un-conserved region are variable. So it can't clarify the hox members base on the ortholog method.
1.2 tblastx
   Firstly, we take the protein and annotation of hox family memebers from the fruit fly genome. We utilizing the highly conserved region search the corresponding hox candiates in Locusta migratoria. And then the Schistocerca species genome and annotation were used to check the hox candates again by tblastx.
   The results are in https://docs.google.com/spreadsheets/d/1kY4OBmMStM12mcQGZ1Ir3BdIPqlKyHbl8GAv9-C__vQ/edit?usp=sharing

  1.2.1 Download the sequences
  wget the Schistocerca species genome from PRJNA772266 and fruit fly Drosophila melanogaster.
  grep the hox genes name from gff files
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
          
  1.2.2 tblastx
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

1.2.3 plot the syteny
1.2.3.1 Plotsr  https://github.com/schneebergerlab/plotsr?tab=readme-ov-file
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

1.2.3.2 SVbyEye https://htmlpreview.github.io/?https://github.com/daewoooo/SVbyEye/blob/master/man/doc/SVbyEye.html#generate-all-versus-all-minimap-alignments
      Easy to use, from minimap2 to plot.
      /home/shangao/software/cactus-bin-v2.6.7/bin/minimap2 -x asm20 -c -eqx -secondary=no fruitfly_NT_033777.3.fasta iqSchSeri2.2_NC_064649.1.fasta > fruitfly2iqSchSeri2.2.align
      R
      library(SVbyEye)
      paf.table <- readPaf(paf.file = paf.file,include.paf.tags = TRUE, restrict.paf.tags = "cg")
      filterPaf(paf.table = paf.table, min.align.len = 100000)
      plotMiro(paf.table = paf.table, color.by = "identity")
      dev.off()
  
        

      
