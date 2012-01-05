'''
This program is actually a variant of `verification_for_bowtie` function  for paired-end data.
'''


import sys
import subprocess
from selected_the_multiply_mapped_reads_of_the_same_genomic_position import *

def pairize(bundle):
	key=bundle.keys()[0]
	bundle=bundle[key]
	if len(bundle)%2==1:
		print >>sys.stderr,bundle
		return False
	
	pairs=[]
	for i in xrange(1,len(bundle),2):
		tab_pos=bundle[i-1].find('\t')
		if bundle[i-1][:tab_pos-2]==bundle[i][:tab_pos-2]:
			pairs.append((bundle[i-1],bundle[i]))
		else:
			print >>sys.stderr,pairs
			return False
	
	return pairs

def translate_into_genomic(records,file1,file2,ENST_structures):
	'''
	file1 is a file handler for outputing result
	file2 is a file handler for outputing incorrect result.
	Calling this function must give both of two file handlers.
	Use /dev/null if one do not want output anything.
	'''
	if len(records)==0:
		return dict()
	
	pairs=pairize(records)
	if pairs==False:
		print >>sys.stderr,"ERROR:pairize(bundle)",records
		return False
	
	pair_info=set()
	for pair in pairs:
		end1=pair[0].split('\t')
		end2=pair[1].split('\t')
		tmp_read_length=adjust_length_by_CIGAR(end1[5],len(end1[9]))
		target1=genome_GPS(end1[2],int(end1[3]),ENST_structures,tmp_read_length)
		if target1==False or target1[1]==False or target1[2]==False:
			print >>sys.stderr,end1
			continue
		tmp_read_length=adjust_length_by_CIGAR(end2[5],len(end2[9]))
		target2=genome_GPS(end2[2],int(end2[3]),ENST_structures,tmp_read_length)
		if target2==False or target2[1]==False or target2[2]==False:
			print >>sys.stderr,end2
			continue
		pair_info.add((tuple(target1),tuple(target2)))
	
	if len(pair_info)==1:
		# All alignments has the same genomic position. 
		pair_info=list(pair_info)
		end1=pairs[0][0].split("\t")
		end2=pairs[0][1].split("\t")
		end1=adjust_alignment_for_strands(ENST_structures,end1,pair_info[0][0])
		end2=adjust_alignment_for_strands(ENST_structures,end2,pair_info[0][1])
		file1.write(end1+'\n')
		file1.write(end2+'\n')
	if len(pair_info)>1:
		print >>file2,pairs[0][0][:-2],pair_info
	
	records=None
	return dict()



def verify_for_paired_end_reads(input_file,output_file,ENST_structures):
	fhr=subprocess.Popen(input_file, shell=True,stdout=subprocess.PIPE).stdout
	fhw=open(output_file,"w+")
	fhw_mul=open(output_file+".mul","w+")
	bundle=dict()
	for line in fhr:
		f=line.split("\t")
		if (int(f[1])& 4)==4:
			continue
		
		ID=f[0][:-2]
		if len(bundle)==0:
			bundle[ID]=[]
		elif ID not in bundle:
			bundle=translate_into_genomic(bundle,fhw,fhw_mul,ENST_structures)
			if bundle==False:
				print >>sys.stderr,bundle,ID,line
			bundle[ID]=[]
		bundle[ID].append(line.strip("\r\n"))
	
	bundle=translate_into_genomic(bundle,fhw,fhw_mul,ENST_structures)
	
	fhw.close()
	fhw_mul.close()
	return True

