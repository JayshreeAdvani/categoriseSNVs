#!/usr/bin/env python
import re,sys

cds_coords = {}
for  line in open(sys.argv[1]):
	line = line.replace("\n","").strip()
	if line.strip() != "":
		line_arr = line.split("\t")
		if len(line_arr) >= 4 and line_arr[2].strip() == "CDS":
			chr = line_arr[0]
			Name = line_arr[8].split("Parent=")[1].split(";")[0].replace('"','').strip()
			for x in range(int(line_arr[3]), int(line_arr[4])+1):
				coord = chr.strip()+"#"+str(x)
				if cds_coords.has_key(coord):
					if Name.strip() not in cds_coords[coord]:
						cds_coords[coord].append(Name.strip())
				else:
					cds_coords[coord]=[Name]
nucleotides = ['A','T','G','C']
for line in open(sys.argv[2]):
	line = line.replace("\n","").strip()
	if line.strip() != "":
		line_arr = line.split("\t")
		if len(line_arr) >= 5:
			variant_allele = line_arr[4]
			if variant_allele in nucleotides:
				chromosome = line_arr[0]
				coord = line_arr[1]
				comb_str = chromosome+"#"+coord
				if cds_coords.has_key(comb_str):
					acc = ";".join(cds_coords[comb_str])
					
					print("%s\t%s\n"%(acc, line.strip()))

