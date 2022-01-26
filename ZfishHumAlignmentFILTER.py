#settings
import os
os.chdir("/project/ktroyer_545/Dana/working_folder/")
#read .bed files using pandas
#module load python
import module load python
import pandas as pd
filename1 = "GSM1817171_human1_K4me1_rep1.bed"
zfish = pd.read_csv(filename1, sep = '\t', header = None, index_col = False)
#visualize the imported dataframes to ensure validity 
#zfish

import subprocess
#Remove promoters and exons
for filename in ./*.narrowPeak; do
	bedtools intersect -v -a "${filename}" -b DanRer.GRCz11.98.promoter.all.bed > "${filename%%.narrowPeak}".NP.narrowPeak
done

for filename in ./*.NP.narrowPeak; do
	bedtools intersect -v -a "${filename}" -b DanRer.GRCz11.98.exon.all.bed > "${filename%%.NP.narrowPeak}".distal.narrowPeak
done

#Filter for peaks within 200kb of HC/SC TSS
sort-bed DanRer.GRCz11.98.TSS.HCvPG.bed > DanRer.GRCz11.98.TSS.HCvPG.sorted.bed
sort-bed DanRer.GRCz11.98.TSS.SCvPG.bed > DanRer.GRCz11.98.TSS.SCvPG.sorted.bed

closest-features --dist --delim '\t' --closest NCC_36hpf_HC_peaks.distal.narrowPeak DanRer.GRCz11.98.TSS.HCvPG.sorted.bed > NCC_36hpf_HC_peaks.distal.HCTSS.narrowPeak
closest-features --dist --delim '\t' --closest NCC_36hpf_SC_peaks.distal.narrowPeak DanRer.GRCz11.98.TSS.HCvPG.sorted.bed > NCC_36hpf_SC_peaks.distal.HCTSS.narrowPeak
closest-features --dist --delim '\t' --closest NCC_36hpf_SC_peaks.distal.narrowPeak DanRer.GRCz11.98.TSS.SCvPG.sorted.bed > NCC_36hpf_SC_peaks.distal.SCTSS.narrowPeak

setwd("G:/My Drive/PhD research/NGS_Data/Mouse ATAC&C&R/PhyloP/Fish_Peaks/")

# read in PhyloP scored peaks with assigned HC/SC TSS
HC.HCTSS.bed <- read.table("NCC_36hpf_HC_peaks.distal.HCTSS.narrowPeak", header = F, sep = "\t", quote = "\"", dec = ".", fill = T)
SC.HCTSS.bed <- read.table("NCC_36hpf_SC_peaks.distal.HCTSS.narrowPeak", header = F, sep = "\t", quote = "\"", dec = ".", fill = T)
SC.SCTSS.bed <- read.table("NCC_36hpf_SC_peaks.distal.SCTSS.narrowPeak", header = F, sep = "\t", quote = "\"", dec = ".", fill = T)

# subset for HC with HCTSS within 200kb
HC.HCTSS.200kb.bed <- subset(HC.HCTSS.bed, subset = abs(HC.HCTSS.bed$V17) < 200000)
# subset for SC with HCTSS within 200kb
SC.HCTSS.200kb.bed <- subset(SC.HCTSS.bed, subset = abs(SC.HCTSS.bed$V17) < 200000)
# subset for SC with SCTSS within 200kb
SC.SCTSS.200kb.bed <- subset(SC.SCTSS.bed, subset = abs(SC.SCTSS.bed$V17) < 200000)


write.table(HC.HCTSS.200kb.bed, file = "NCC_36hpf_HC_peaks.distal.HCTSS200kb.narrowPeak", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(SC.HCTSS.200kb.bed, file = "NCC_36hpf_SC_peaks.distal.HCTSS200kb.narrowPeak", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(SC.SCTSS.200kb.bed, file = "NCC_36hpf_SC_peaks.distal.SCTSS200kb.narrowPeak", sep = "\t", row.names = F, col.names = F, quote = F)

# Save a version of PhyloP scored peaks within 200kb of HC/SC TSS, without the genome coordinates of the TSS for visualization in genome browser
library(dplyr)
HC.HCTSS.200kb.genenameonly.bed <- HC.HCTSS.200kb.bed %>% select(V1:V10, V14, V17)
SC.HCTSS.200kb.genenameonly.bed <- SC.HCTSS.200kb.bed %>% select(V1:V10, V14, V17)
SC.SCTSS.200kb.genenameonly.bed <- SC.SCTSS.200kb.bed %>% select(V1:V10, V14, V17)

write.table(HC.HCTSS.200kb.genenameonly.bed, file = "NCC_36hpf_HC_peaks.distal.HCTSS200kb.genenameonly.narrowPeak", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(SC.HCTSS.200kb.genenameonly.bed, file = "NCC_36hpf_SC_peaks.distal.HCTSS200kb.genenameonly.narrowPeak", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(SC.SCTSS.200kb.genenameonly.bed, file = "NCC_36hpf_SC_peaks.distal.SCTSS200kb.genenameonly.narrowPeak", sep = "\t", row.names = F, col.names = F, quote = F)
