#!/usr/bin/env python
import re,sys
f1 = open("GCF_000195955.2_ASM19595v2_genomic.fasta","r")
t1 = f1.read()
f1.close()
t1_arr = t1.split(">")
genome_seq = {}
for x in t1_arr[1:]:
	xarr = x.split("\n",1)
	acc = xarr[0].split(" ")[0].split("|")[3]
	print acc
	seq = xarr[1].replace("\n","").strip()
	genome_seq[acc]=seq.strip()
		
def fetch_SNV(snp, coords, strand,chr,change_nuc,acc):
	n ={"A":"T","T":"A","C":"G","G":"C","N":"N"}
	codon_dict_1={'TTT':'F','TCT':'S','TAT':'Y','TGT':'C','TTC':'F','TCC':'S','TAC':'Y','TGC':'C','TTA':'L','TCA':'S','TAA':'*','TGA':'*','TTG':'L','TCG':'S','TAG':'*','TGG':'W','CTT':'L','CCT':'P','CAT':'H','CGT':'R','CTC':'L','CCC':'P','CAC':'H','CGC':'R','CTA':'L','CCA':'P','CAA':'Q','CGA':'R','CTG':'L','CCG':'P','CAG':'Q','CGG':'R','ATT':'I','ACT':'T','AAT':'N','AGT':'S','ATC':'I','ACC':'T','AAC':'N','AGC':'S','ATA':'I','ACA':'T','AAA':'K','AGA':'R','ATG':'M','ACG':'T','AAG':'K','AGG':'R','GTT':'V','GCT':'A','GAT':'D','GGT':'G','GTC':'V','GCC':'A','GAC':'D','GGC':'G','GTA':'V','GCA':'A','GAA':'E','GGA':'G','GTG':'V','GCG':'A','GAG':'E','GGG':'G'}

	all_coords = []
	snp= int(snp) 
	for x in coords:
		c = x.split("-")
		if strand.strip() == "-":
			c1 = range(int(c[1]), int(c[0])-1,-1)
		else:
			c1 = range(int(c[0]),int(c[1])+1)
		all_coords.extend(c1)
	seq_pos = all_coords.index(int(snp))+1
	aa_pos = seq_pos/3
	remainder = seq_pos%3
	seq_pos = seq_pos -1
	g_seq = genome_seq[chr.strip()].upper().strip()
	if remainder==0:
		codon = [g_seq[all_coords[seq_pos-2]-1],g_seq[all_coords[seq_pos-1]-1],g_seq[all_coords[seq_pos]-1]]
		codon2 = [g_seq[all_coords[seq_pos-2]-1],g_seq[all_coords[seq_pos-1]-1],change_nuc]
	if remainder == 1:
		codon = [g_seq[all_coords[seq_pos]-1],g_seq[all_coords[seq_pos+1]-1],g_seq[all_coords[seq_pos+2]-1]]
		codon2 = [change_nuc,g_seq[all_coords[seq_pos+1]-1],g_seq[all_coords[seq_pos+2]-1]]
		aa_pos = aa_pos +1
	if remainder == 2:
		codon = [g_seq[all_coords[seq_pos-1]-1],g_seq[all_coords[seq_pos]-1],g_seq[all_coords[seq_pos+1]-1]]
		codon2 = [g_seq[all_coords[seq_pos-1]-1],change_nuc,g_seq[all_coords[seq_pos+1]-1]]
		aa_pos = aa_pos +1
	if strand == "-":
		codon[0]=n[codon[0]]
		codon[1]=n[codon[1]]
		codon[2]=n[codon[2]]
		codon2[0]=n[codon2[0]]
		codon2[1]=n[codon2[1]]
		codon2[2]=n[codon2[2]]
	codon = "".join(codon)
	codon2 = "".join(codon2)
	if codon_dict_1.has_key(codon) and codon_dict_1.has_key(codon2):
		aa1 = codon_dict_1[codon]
		aa2 = codon_dict_1[codon2]
		if aa1.strip() != aa2.strip():
			return_str = chr+":"+strand+":"+g_seq[all_coords[seq_pos]-1]+str(snp)+change_nuc+":"+codon+"->"+codon2+":"+acc+":"+aa1+str(aa_pos)+aa2
			print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(chr, strand, g_seq[all_coords[seq_pos]-1],str(snp),change_nuc,codon,codon2,acc,aa1, str(aa_pos),aa2,return_str))
		else:
			return_str=''
	else:
		return_str=''
	return return_str
		
transcript_coords = {}
transcript_strand = {}
transcript_chr = {}
for  line in open("Mycobacterium_tuberculosis_H37Rv.gtf"):
	line = line.replace("\n","").strip()
	if line.strip() != "":
		line_arr = line.split("\t")
		if len(line_arr) >= 4 and line_arr[2].strip() == "CDS":
			chr = line_arr[0]
			transcript_id = line_arr[8].split("transcript_id")[1].split(";")[0].replace('"','').strip()
			coords = line_arr[3]+"-"+line_arr[4]
			transcript_strand[transcript_id]= line_arr[6]
			transcript_chr[transcript_id]=line_arr[0].strip()
			if transcript_coords.has_key(transcript_id):
				transcript_coords[transcript_id].append(coords)
			else:
				transcript_coords[transcript_id]=[coords]

nuc =["A","T","C","G"]
for line in open(sys.argv[1]):
	line = line.replace("\n","").strip()
	if line.strip() != "":
		line_arr = line.split("\t")
		m_acc = line_arr[0].strip().split(";")
		for acc in m_acc:
			coords = transcript_coords[acc.strip()]
			strand = transcript_strand[acc.strip()]
			chr = transcript_chr[acc.strip()]
			snp = line_arr[2]
			change_nuc = line_arr[5].strip().upper()
			if change_nuc in nuc:
				snp_codon = fetch_SNV(snp,coords,strand,chr,change_nuc,acc)
