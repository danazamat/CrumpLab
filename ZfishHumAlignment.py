##settings
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


#iterate over rows in zfish enhancer file to 
#identify alignment in mouse genome near same gene
import subprocess
i = 0
for row in zfish.itertuples():
	#generate coordinates for enhancers in each row and store in temp file
	zfish_chr, zfish_start, zfish_end = row[1:4]
	zfish_coord = str(zfish_chr) + "\t" + str(zfish_start) + "\t" + str(zfish_end) + "\n"
	bed = open("zfish_coord.bed", "w+")
	bed.write(zfish_coord)
	bed.close()
	#zfish_chr_fa = "./mm_chr/" + str(zfish_chr) + ".fa"

	# fasta from zfish
	zfish_seq = subprocess.run(['bedtools', 'getfasta', '-fi', '/project/ktroyer_545/Dana/ref_genome/danRer11.fa', '-bed', 'zfish_coord.bed'], stdout = subprocess.PIPE, universal_newlines = True)
	fasta = open("zfish.fa", "w+")
	fasta.write(zfish_seq.stdout)
	fasta.close()

	# perform local alignment
	# return alignments with evalue less than 0.01
	# sseqid sstart send are coordinates for mouse sequence 
	align = subprocess.run(['blastn', '-task', 'blastn', '-query', 'zfish.fa', '-subject', 'mm10.fa' , '-evalue', '0.01', '-outfmt', '6 qseqid sseqid sstart send score length pident evalue'], stdout = subprocess.PIPE, universal_newlines = True)

	# store alignment info into new text --> appending w each loop
	if len(align.stdout) > 0:
		alignstdout = open("align.stdout.txt", "w+")
		alignstdout.write(align.stdout)
		alignstdout.close()
		blastn = pd.read_csv("align.stdout.txt", sep = '\t', header = None, index_col = False)

		# check what nearest gene for alignment 
		for row2 in blastn.itertuples():
			mouse_chr, mouse_start, mouse_end = row2[2:5]
			mouse_coord = str(mouse_chr) + '\t' + str(mouse_start) +'\t' + str(mouse_end) + '\n'
			bed2 = open("mouse_coord.bed", "w+")
			bed2.write(mouse_coord)
			bed2.close()

			# use bedops to locate nearest to aligned sequence
			# TSS file has all mouse gene names to be converted to zfish orthologous gene names for comparison 
			bedops = subprocess.run(['closest-features', '--dist', '--delim', '\t', '--closest', 'mouse_coord.bed', 'Mus_musculus.GRCm38.101.TSS.HCvPG.sorted.bed'], 
				stdout = subprocess.PIPE, universal_newlines = True)
			closest_gene = open("closest_gene.bed", "w+")
			closest_gene.write(bedops.stdout)
			closest_gene.close()
			closestgene = pd.read_csv("closest_gene.bed", sep = '\t', header = None, index_col = False)

			#bedops will return NA if no genes near 
			if pd.isna(closestgene.iloc[0,4]) == True:
				continue

			# if bedops assigned gene is the same as zfish peak and within 200kb, write output
			if (closestgene.loc[0,6] == row[11]) & (abs(closestgene.loc[0,9]) <= 200000):
				zfish_coord = row2[1]
				score, length, pident, evalue = row2[5:9]
				output = str(row[11]) + '\t' + str(zfish_coord) + '\t' + str(mouse_chr) + '\t' + str(mouse_start) + '\t' + str(mouse_end) + '\t' + str(score)
				 + '\t' + str(length) + '\t' + str(pident) + '\t' + str(evalue) + '\n' 
				alignment = open("mmAlignmentZfish.txt", "a+")
				alignment.write(output)
				alignment.close()
	i = i + 1
	print(str(i) + " check!")      


#GENERAL FILTER: 
#Remove promoters and exons
for filename in ./*.narrowPeak; do
	bedtools intersect -v -a "${filename}" -b Mus_musculus.GRCm38.101.TSS.HCvPG.sorted.bed > "${filename%%.narrowPeak}".NP.narrowPeak

for filename in ./*.NP.narrowPeak; do
	bedtools intersect -v -a "${filename}" -b Mus_musculus.GRCm38.101.TSS.HCvPG.sorted.bed > "${filename%%.NP.narrowPeak}".distal.narrowPeak

#Filter for peaks within 200kb
sort-bed Mus_musculus.GRCm38.101.TSS.HCvPG.sorted.bed > Mus_musculus.GRCm38.101.TSS.HCvPG.sorted.bed
sort-bed Mus_musculus.GRCm38.101.TSS.HCvPG.sorted.bed > Mus_musculus.GRCm38.101.TSS.HCvPG.sorted.bed


#FILTER FOR S LINE: 
# S = sequence of zebrafish (!!) --> align from this file to. 
# If there is a match then return the S line with the human coordinates. return the human sequence.
#can return every S line and then trim the file later on for human coordinates 






#final output should be human coordinates****--> human in one line and zebrafish on another 
# SEPARATE 





