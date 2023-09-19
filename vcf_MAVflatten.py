#!/usr/bin/env python
import sys
import re
from utils import return_file_handle

"""
The script separates/flattens MAVs in a given VCF and keeps the Bi-allelic variants intact.
These are the modifications which the script does
	Reg GT		- It splits the Genotype fields GT, AD, DP, GQ, PL
	Reg AC, AN, AF	- Modified the GT field. Recalculates AC, AN, AF. 
	Reg GQ 		- Keeps the same for all ALT of an MAV
	Reg DP 		- Keeps the sum of AD for each ALT allele of an MAV
	Reg AD		- Takes the AD values for each ALT allele
	Reg PL 		- Takes the PL value for each ALT allele

The script doesnot remove AC=0 variants. Use the script /cgrl/CommonDATA/Python_Scripts/vcf_RemoveACzerovariants.py to remove AC=0 variants.
"""

vcf = sys.argv[1]

PL_index = {}

PL_index = {
	"1"  : [0,1,2],
	"2"  : [0,3,5],
	"3"  : [0,6,9],
	"4"  : [0,10,14],
	"5"  : [0,15,20],
	"6"  : [0,21,27],
	"7"  : [0,28,35],
	"8"  : [0,36,44],
	"9"  : [0,45,54],
	"10" : [0,55,65]
	}
fvcf = return_file_handle(vcf)
for line in fvcf:
	line = line.strip().split('\t')
	if re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
		try:
			GT_pos = line[8].split(':').index('GT')
		except:
			GT_pos = ''
		try:
			DP_pos = line[8].split(':').index('DP')
		except:
			DP_pos = ''
		try:
			AD_pos = line[8].split(':').index('AD')
		except:
			AD_pos = ''
		try:
			GQ_pos = line[8].split(':').index('GQ')
		except:
			GQ_pos = ''
		try:
			PL_pos = line[8].split(':').index('PL')
		except:
			PL_pos = ''

		formatField = line[8].split(':')	# extracting the formatField values

		"""
		Looping through each ALT allele
		"""
		if ',' in line[4] :
			final_genotype = []
			
			for alt in range(0,len(line[4].split(','))):
				GT_new = []		
				alt_flag = 0
				new_GT = []
				new_AD = []
				new_genotype = []
				final_genotype = []
				temp_alt = alt + 1 #actual ALT GT values
				if temp_alt == 1:
					pl_start = 0
					pl_stop = 2
				else:
					pl_start = pl_stop + 1
					pl_stop = pl_start + temp_alt
				for genotype_fields in line[9:]:	#manipulating the GT for MAVs
					"""
					Processing of AD fields
					"""
					if 'AD' not in formatField :
						ad = '.'
					elif 'AD' in formatField and len(genotype_fields.split(':')) <= AD_pos :
						ad = '.'
					elif len(genotype_fields.split(':')) > AD_pos and genotype_fields.split(':')[AD_pos] == '.' :
						ad = '.'
					elif len(genotype_fields.split(':')) > AD_pos and genotype_fields.split(':')[AD_pos] != '.' :		#Getting AD and DP values for each ALT allele
						ad1 = genotype_fields.split(':')[AD_pos].split(',')[0]
						ad2 = genotype_fields.split(':')[AD_pos].split(',')[temp_alt]
						ad = ad1 + ',' + ad2
					else:
						ad = '.'
					"""
					Processing of DP fields
					"""
					if 'DP' not in formatField :
						dp = '.'
					elif 'DP' in formatField and len(genotype_fields.split(':')) <= DP_pos :
						dp = '.'
					elif len(genotype_fields.split(':')) > DP_pos and genotype_fields.split(':')[DP_pos] == '.' :
						dp = '.'
					elif len(genotype_fields.split(':')) > DP_pos and genotype_fields.split(':')[DP_pos] != '.' :		#Getting GQ values for each genotype
						dp = genotype_fields.split(':')[DP_pos]
					else:
						dp = '.'
					"""
					Processing of GQ fields
					"""
					if 'GQ' not in formatField :
						gq = '.'
					elif 'GQ' in formatField and len(genotype_fields.split(':')) <= GQ_pos :
						gq = '.'
					elif len(genotype_fields.split(':')) > GQ_pos and genotype_fields.split(':')[GQ_pos] == '.' :
						gq = '.'
					elif len(genotype_fields.split(':')) > GQ_pos and genotype_fields.split(':')[GQ_pos] != '.' :		#Getting GQ values for each genotype
						gq = genotype_fields.split(':')[GQ_pos]
					else:
						gq = 0
					# Mnaipulating GT and PL fields
					"""
					Processing of GT fields
					"""
					temp_gt = []
					val_temp = []
					if 'GT' not in formatField :
						temp_gt = ['.', '.']
					elif 'GT' not in formatField and len(genotype_fields.split(':')) <= GT_pos :
						temp_gt = ['.', '.']
					elif len(genotype_fields.split(':')) > GT_pos and genotype_fields.split(':')[GT_pos] == '.' :
						temp_gt = ['.', '.']
					elif len(genotype_fields.split(':')) > GT_pos  and genotype_fields.split(':')[GT_pos] != '.' :
						if "/" in genotype_fields.split(':')[GT_pos]:
							for val in genotype_fields.split(':')[GT_pos].split('/'):
								alt_flag = 1
								val_temp.append(val)
								if val == str(temp_alt) and str(val) != '0':
									if val == '1':
										temp_gt.append(val)
									else:
										temp_gt.append('1')
								elif val != str(temp_alt) and str(val) != '0':
									new_val = val.replace(val,".")
									temp_gt.append(new_val)
								elif val == '0':
									temp_gt.append(val)
								elif val == ".":
									temp_gt.append(val)
						elif "|" in genotype_fields.split(':')[GT_pos]:
							for val in genotype_fields.split(':')[GT_pos].split('|'):
								alt_flag = 1
								val_temp.append(val)
								if val == str(temp_alt) and str(val) != '0':
									if val == '1':
										temp_gt.append(val)
									else:
										temp_gt.append('1')
								elif val != str(temp_alt) and str(val) != '0':
									new_val = val.replace(val,".")
									temp_gt.append(new_val)
								elif val == '0':
									temp_gt.append(val)
								elif val == ".":
									temp_gt.append(val)
					else:
						temp_gt = ['.','.']
					
					GT_new = GT_new + temp_gt
					"""
					Processing of PL fields
					"""
					if '.' in temp_gt :
						aa = '.'
					else:
						aa =  max(val_temp)
					if 'PL' not in formatField :
						PL_str = '.'
					elif PL_pos != '' and len(genotype_fields.split(':')) <= PL_pos :
						PL_str = '.'
					else:
						if aa == '0':
							aa  = '1'
							PL_1 = PL_index[aa][0]
							PL_2 = PL_index[aa][1]
							PL_3 = PL_index[aa][2]
	                	                        # PL_2 = PL_index[aa][1]
        	                	                # PL_3 = PL_index[aa][2]
							PL_str = genotype_fields.split(':')[PL_pos].split(',')[PL_1] + ',' +  genotype_fields.split(':')[PL_pos].split(',')[PL_2] + ',' +  genotype_fields.split(':')[PL_pos].split(',')[PL_3]
						elif aa == '.':
							PL_str = '.'
						
						else:
							PL_1 = PL_index[aa][0]
							PL_2 = PL_index[aa][1]
							PL_3 = PL_index[aa][2]
	                        	                # PL_2 = PL_index[aa][1]
        	                        	        # PL_3 = PL_index[aa][2]
							PL_str = genotype_fields.split(':')[PL_pos].split(',')[PL_1] + ',' +  genotype_fields.split(':')[PL_pos].split(',')[PL_2] + ',' +  genotype_fields.split(':')[PL_pos].split(',')[PL_3]
					final_genotype_str = "/".join(temp_gt) + ':' + str(ad) + ':' + str(dp) + ':' + str(gq) + ':' + str(PL_str) # the final gentoype field for each sampleID		
					final_genotype.append(final_genotype_str) # appending the genotype fields for each sampleID to a list which is kept for each variant
				temp_an = 0
				an_new = ''
				ac_new = ''
				af_new = ''
				temp_an = int(GT_new.count('0')) + int( GT_new.count('1'))
				an_new = 'AN=' + str(temp_an)
				ac_new = 'AC=' + str(GT_new.count('1'))
				#print(ac_new,an_new)
				try:
					af_new = 'AF=' + str(int(GT_new.count('1'))/float(temp_an))
				except ZeroDivisionError:
					continue
					#print(final_genotype)
					#break
				for inf in line[7].split(';'):
					if inf.startswith('AC='):
						line[7] = line[7].replace(inf,ac_new)
					elif inf.startswith('AN='):
						line[7] = line[7].replace(inf,an_new)
					elif inf.startswith('AF='):
						line[7] = line[7].replace(inf,af_new)
						
				print '\t'.join(line[0:4] + [line[4].split(',')[alt]] + line[5:8] + [":".join(['GT','AD','DP','GQ','PL'])] + final_genotype)
		else:
			print '\t'.join(line)	#if the variant is not a MAV, just printing out the values as it is in the VCF
	else:	# printing the non variant lines in the VCF
		print '\t'.join(line)
fvcf.close()
