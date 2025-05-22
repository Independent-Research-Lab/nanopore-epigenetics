# nanopore-epigenetics
#Human Ageing Epigenetics Using Oxford Nanopore Long Reads DNA RRMS Sequencing

This script uses DNA methylation .bed files obtained with Oxford Nanopore human DNA sequencing to find the CpGs methylated in direct or reversed correlattion with ageing.  

To run this R script, you need to obtain methypation .bed files following the below steps:

1. Get your nanopore data basecalled for each sample, using dorado from MinKnow interface with the model HAC. Don't forget to include modifications 5mc and 5hmc! You also need a linux machine or WSL running any linux on which you must install samtools and modkit, through yum or apt or your linux package manager or from the below repositories:

https://github.com/samtools/
https://github.com/nanoporetech/modkit

3. The resulted *.bam files, in our case located in the folder basecalling/pass/  ,must be concatenated using samtools, like this:

samtools merge merged.bam basecalling/pass/*.bam

3. Do some check, sorting and indexing of the resulted .bam file like this:

samtools quickcheck merged.bam; samtools sort -o merged.sorted.bam merged.bam; samtools index merged.sorted.bam

4. do pileup with modkit:

modkit pileup merged.sorted.bam pileup_mod.bed

5. remane your bam file in a special format which our script recognizes as input. That format must be something like xxNN, where xx are any 2 alphabetical non numerical characters, and NN are 2 numerical digits representing the age of the dample donor. If the sample is from a 76 years old person, rename that as xx76.bam:
mv merged.sorted.bam xx76.bam

6. do points 1-5 for all your samples

7. Copy or move all the resulting *.bam files in a folder and set that as R working directory.
for convenience and reproducibility we included our .bam files in the link below.
8. Due to limitations on file size, our data can not be uploaded on github. It is available on the link below instead, as complete content of R working directory: 

http://gofile.me/76Mv5/imaOtyFz3


#To support our research work, feel free to donate any amount here:

https://www.independent-research.ro/pages/doneaza-cu-cardul

#Good luck with your epigenetics and ageing research!
