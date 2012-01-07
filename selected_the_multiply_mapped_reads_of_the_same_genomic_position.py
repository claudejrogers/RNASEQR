'''
This program is to identify those multiply mapped reads which are mapped to different isoform of the same genes.

This program has three steps.
1. Load the realtions of ENSTs and ENSGs.
2. Select the reads which maps to different ENSTs of the same ENSG
3. Dump them to the output file

The Ensembl Annotation lists exons from 5' to 3'. The start and end positions of each exon is based on coordinates of positive strand. If the transcript is on the reverse strand, then the start and end positions are based on the coordinates of positive strand. But the DNA sequence in the cDNA file is from the reverse strand.

Note that the reported position of alignment in SAM format is always 1-based coordiantes of refernce strand. (Although we will assume the reference is postive strand usually. It is not always ture when the reference is transcriptome.) Not like BLAST/BLAT, they will report algined positions by the head of query string. If the reversed query string aligned to the reference string, then BLAST will report the position of end of query string.

E.g., the query string is CTGAAGAAGT, and its complementary string aligns to the reference.
CGGCAGTAGGCAATAAACTTCTTCAGCTTGGCCAGGTCAATCTCGCCCTCCGCAGC
                ACTTCTTCAG
The Blast will report 26, and SAM will be 17.

If the gene is on the reversed strand, then it become more complex.

CGGCAGTAGGCAATAAACTTCTTCAGCTTGGCCAGGTCAATCTCGCCCTCCGCAGC
                ACTTCTTCAG
The Blast will report 31, and SAM will be 40.

'''

import re
from optparse import OptionParser
from array import *
from interval_tree import *
from map_and_tables import *
from sys import stdout
import subprocess

chr_boundary=dict(zip(chr_id,chr_length))
chr_id_name_map=dict(zip(chr_id,chr_name))
chr_name_id_map=dict(zip(chr_name,chr_id))


SAM_XA_parser = re.compile(r"([^:;,]+,[\+|\-][0-9]+,[^;]+,[^;]+);")
ENSG_parser = re.compile('gene_id \"([^\"]*)\";')
ENST_parser = re.compile('transcript_id \"([^\"]*)\";')
exon_number_parser = re.compile('exon_number \"([0-9]*)\";')
CIGAR_parser = re.compile(r'([0-9]+[M|I|D|N|S|H|P])')
Number_parser = re.compile(r'([0-9]+)')
Insertion_parser = re.compile(r'([0-9]+)I')
Deletion_parser = re.compile(r'([0-9]+)D')
Skip_parser = re.compile(r'([0-9]+)N')
SoftClip_parser = re.compile(r'([0-9]+)S')

def Is_SAM_header(line,extension=[]):
	SAM_spec=["@HD","@SQ","@RG","@PG","@CO"]+extension
	SAM_spec=set(SAM_spec)
	
	if line[:3] in SAM_spec:
		return True
	return False 

def extract_SAM_XA(line):
	XA_pos=line.find("XA:Z:")
	if XA_pos == -1:
		return []
	result=line[XA_pos+5:].strip(" \r\n\t").split(";")
	j=len(result)
	for i in xrange(len(result)):
		'''
		if 3!=result[i].count(","):
			j=i
			break
		'''
		if  (result[i]=='' or (-1!=result[i].find(" ")) or (-1!=result[i].find("\t"))):
			j=i
			break
	return result[:j]


def get_exon_width(exon_number,ENST):
	# ENST is ENST_structures[ENST_id]
	size=len(ENST[1])
	if exon_number>=size:
		return False
	elif exon_number==0:
		return ENST[1][0]
	else:
		return ENST[1][exon_number]-ENST[1][exon_number-1]


def extract_exon_info(line):
	global ENSG_parser
	global ENST_parser
	global exon_number_parser
	a=ENSG_parser.search(line)
	b=ENST_parser.search(line)
	c=exon_number_parser.search(line)
	if a!=None:
		a=a.group(1)
	if b!=None:
		b=b.group(1)
	if c!=None:
		c=c.group(1)
	
	return (a,b,c)

def extend_ENST_info_records(gtf_file,ENST_structures):
	global ENST_parser
	f=open(gtf_file,"r")
	for line in f:
		li=line.split("\t")
		if li[1]!="protein_coding":
			continue
		tmp=ENST_parser.findall(line)
		if tmp == [] or len(tmp)>1:
			return False
		else:
			tmp=tmp[0]
			if tmp in ENST_structures:
				ENST_structures[tmp]+=("C",)
	f.close()
	for ENST in  ENST_structures.iterkeys():
		if len(ENST_structures[ENST])<5:
			ENST_structures[ENST]+=("NC",)
	return True

def extract_gene_isoform_relation(gtf_file):
	ENST2ENSG=dict()
	fhr=open(gtf_file,"r")
	for line in fhr:
		line_tmp=line.strip(" \r\n\t")
		line_tmp=line_tmp.split("\t")
		# we only consider exons
		if line_tmp[2]=="exon":
			tmp_ENSG,tmp_ENST,tmp_exon_number=extract_exon_info(line)
			if tmp_ENST not in ENST2ENSG:
				ENST2ENSG[tmp_ENST]=tmp_ENSG
		
	fhr.close()
	return ENST2ENSG

def complement_sequence(seq):
	table = 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'\
                'xTxGxxxCxxxxxxNxxxxxAxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'\
                'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'\
                'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
	comp = seq[::-1].translate(table)
	if 'x' in comp:
		raise ValueError, 'Wrong Sequence Format!'
	return comp

def reverse_CIGAR(CIGAR):
	CIGAR=CIGAR_parser.findall(CIGAR)
	CIGAR.reverse()
	return "".join(CIGAR)


def extract_gene_structure_from_bed(annotation_file):
	'''
	BED format
	http://genome.ucsc.edu/FAQ/FAQformat#format1
	This format use zero-based semi-open notaton, i.e., [a,b).
	It means the b-th base is not included in the interval.
	
	The BED format is not defined very precisely. Thus it is not
	define the order of exons for each isoform. Hence, this program
	is based on the output of UCSC GB.
	http://genome.ucsc.edu/cgi-bin/hgTables 	(choose output format as BED)
	'''
	global chr_id
	chr_id_tmp=[]
	for i in chr_id:
		if i=="MT":
			chr_id_tmp.append("chrM")
		else:
			chr_id_tmp.append("chr"+i)
	mapping=dict(zip(chr_id_tmp,chr_id))
	
	ENST_structures=dict()
	fhr=open(annotation_file,"r")
	for line in fhr:
		line_tmp=line.strip(" \r\n\t")
		line_tmp=line_tmp.split("\t")
		chrom=line_tmp[0]
		chromStart=int(line_tmp[1])
		if chrom not in mapping:
			continue
		
		try:
			ref_name=line_tmp[3]
			strand=line_tmp[5]
			exon_number=int(line_tmp[9])
			exon_sizes=line_tmp[10].split(",")
			exon_starts=line_tmp[11].split(",")
		except IndexError:
			print "IndexError A:",line
			continue

		if len(exon_sizes)==len(exon_starts):
			for i in xrange(len(exon_sizes)):
				if exon_sizes[i]=='':
					del exon_sizes[i]
				else:
					exon_sizes[i]=int(exon_sizes[i])
				if exon_starts[i]=='':
					del exon_starts[i]
				else:
					exon_starts[i]=int(exon_starts[i])+chromStart+1
			if strand=="-":
				exon_sizes.reverse()
				exon_starts.reverse()
		else:
			print "exon_sizes and exon_starts do not agreee"
			continue
		
		if len(exon_sizes)==exon_number:
			if ref_name not in ENST_structures:
				for i in xrange(1,exon_number):
					exon_sizes[i]=exon_sizes[i]+exon_sizes[i-1]
				ENST_structures[ref_name]=(array("L",exon_starts),array("L",exon_sizes),strand,mapping[chrom])
			else:
				print "Duplicate Transcript", line
				continue
		else:
			print "exon numbers disagrees"
			continue
	return ENST_structures

def extract_gene_structure_from_gtf(annotation_file):
	'''
	GTF2.2: A Gene Annotation Format
	http://mblab.wustl.edu/GTF22.html
	'''
	ENST_structures=dict()
	prev_ENST=None
	prev_ENST_strand=None
	prev_ENST_chromosme=None
	curr_block=dict()
	
	fhr=open(annotation_file,"r")
	for line in fhr:
		line_tmp=line.strip(" \r\n\t")
		line_tmp=line_tmp.split("\t")
		
		# we only consider exons
		if line_tmp[2]=="exon":
			# we only consider 25 chromsomes
			if line_tmp[0] not in chr_name_id_map:
				continue
			
			tmp_ENSG,tmp_ENST,tmp_exon_number=extract_exon_info(line)
			tmp_exon_number=int(tmp_exon_number)
			
			if prev_ENST==tmp_ENST:
				prev_ENST_strand=line_tmp[6]
				prev_ENST_chromosme=line_tmp[0]
				if tmp_exon_number not in curr_block:
					curr_block[tmp_exon_number]=(int(line_tmp[3]),int(line_tmp[4]))
				else:
					print prev_ENST
					print "Two the same exons !!!", line
					print "tmp_ENST,tmp_ENSG,tmp_exon_number",tmp_ENST,tmp_ENSG,tmp_exon_number
					print curr_block
					exit(2)
			else:
				'''
				Enter new ENST block of some gene. 
				'''
				try:
					#Dump the old one first.
					tmp_len=len(curr_block)
					exon_widths=array("L")
					exon_starts=array("L")
					if tmp_len>0:
						for w in xrange(1,tmp_len+1):
							exon_widths.append(curr_block[w][1]-curr_block[w][0]+1)
							exon_starts.append(curr_block[w][0])
						
						for w in xrange(1,tmp_len):
							exon_widths[w]=exon_widths[w]+exon_widths[w-1]
						
						if prev_ENST not in ENST_structures:
							ENST_structures[prev_ENST]=(exon_starts,exon_widths,prev_ENST_strand,prev_ENST_chromosme)
						else:
							print "Two the same ENST !!!",line
							exit(2)
					
					# Initialization for new block
					prev_ENST=tmp_ENST
					prev_ENST_strand=line_tmp[6]
					prev_ENST_chromosme=line_tmp[0]
					curr_block=None
					curr_block=dict()
					curr_block[tmp_exon_number]=(int(line_tmp[3]),int(line_tmp[4]))
				except KeyError:
					print "KeyError:",line
					print "tmp_ENST,tmp_ENSG,tmp_exon_number",tmp_ENST,tmp_ENSG,tmp_exon_number
					print curr_block
					exit(2)
	fhr.close()
	
	# If the last block does not dump yet, dump it!
	tmp_len=len(curr_block)
	exon_widths=array("L")
	exon_starts=array("L")
	if tmp_len>0:
		for w in xrange(1,tmp_len+1):
			exon_widths.append(curr_block[w][1]-curr_block[w][0]+1)
			exon_starts.append(curr_block[w][0])
		
		for w in xrange(1,tmp_len):
			exon_widths[w]=exon_widths[w]+exon_widths[w-1]
		
		if prev_ENST not in ENST_structures:
			ENST_structures[prev_ENST]=(exon_starts,exon_widths,prev_ENST_strand,prev_ENST_chromosme)
		else:
			print "Two the same ENST !!!",line
			exit(2)
		
	return ENST_structures


def extract_gene_structure(annotation_file,file_type="GTF"):
	if file_type=="GTF":
		return extract_gene_structure_from_gtf(annotation_file)
	elif file_type=="BED":
		return extract_gene_structure_from_bed(annotation_file)
	else:
		return False


def compact_CIGARs(old_CIGAR):
	'''
	This function is to compact the CIGAR strings, e.g.,
	it converts 3M30M500D20D7M1000D30M5M to 33M529D7M1000D35M.
	'''
	Cparser = re.compile(r'([+|-]*[0-9]+[M|I|D|N|S|H|P])')
	old_CIGAR=Cparser.findall(old_CIGAR)
	for i in xrange(1,len(old_CIGAR)):
		first,second=old_CIGAR[i-1:i+1]
		if first[-1]==second[-1]:
			old_CIGAR[i]=str(int(first[:-1])+int(second[:-1]))+second[-1]
			old_CIGAR[i-1]=""
	return "".join(old_CIGAR)


def merge_CIGARs(transcript_CIGAR,genome_CIGAR):
	'''
	If we align reads to transcriptome, then we will have transcript_CIGAR, however this aligned region will have a genome_GICAR according to gene structure. Thus we need to merge these two CIGARs.
	'''
	
	global CIGAR_parser,Number_parser
	# There is no given CIGAR string, thus we just return the CIGAR.
	t_structure=CIGAR_parser.findall(transcript_CIGAR)
	g_structure=CIGAR_parser.findall(genome_CIGAR)
	#m_structure=[]
	output_CIGAR=[]
	pop_g_structure=True
	pop_t_structure=True
	g_type='X'
	t_type='X'
	tmp=0
	try:
		while(not (g_type=='' and t_type=='')):
			if (pop_g_structure==True):
				if len(g_structure)>0:
					i=g_structure.pop(0)
					g_type=i[-1:]
					g_value=int(i[:-1])
					pop_g_structure=False
				else:
					g_type=''
			
			if (pop_t_structure==True):
				if len(t_structure)>0:
					j=t_structure.pop(0)
					t_type=j[-1:]
					t_value=int(j[:-1])
					pop_t_structure=False
				else:
					t_type=''
			
			'''
			This condition is only for that transcript_CIGAR is generated by aligner non supporting long gap alignment.
			'''
			if t_type=='M' and g_type=='M':
				tmp=min(g_value,t_value)
				output_CIGAR.append(str(tmp)+'M')
				g_value-=tmp
				t_value-=tmp
			elif t_type=='I' and g_type=='M':
				output_CIGAR.append(str(t_value)+'I')
				t_value=0
			elif t_type=='D' and g_type=='M':
				tmp=min(g_value,t_value)
				output_CIGAR.append(str(tmp)+'D')
				t_value-=tmp
				g_value-=tmp
			elif t_type=='N' and g_type=='M':
				tmp=min(g_value,t_value)
				output_CIGAR.append(str(tmp)+'N')
				t_value-=tmp
				g_value-=tmp
			elif t_type=='S' and g_type=='M':
				output_CIGAR.append(str(t_value)+'S')
				t_value=0
				g_value-=t_value
			elif t_type=='M' and g_type=='N':
				output_CIGAR.append(str(g_value)+'N')
				g_value=0
			elif t_type=='I' and g_type=='N':
				output_CIGAR.append(str(t_value)+'I')
				output_CIGAR.append(str(g_value)+'N')
				g_value=0
				t_value=0
			elif t_type=='D' and g_type=='N':
				output_CIGAR.append(str(g_value)+'N')
				g_value=0
			elif t_type=='N' and g_type=='N':
				output_CIGAR.append(str(g_value+t_value)+'N')
				g_value=0
				t_value=0
			elif t_type=='S' and g_type=='N':
				output_CIGAR.append(str(t_value)+'S')
				t_value=0
				g_value-=g_value
			elif t_type=='' and g_type=='':
				pass
			elif t_type=='':
				output_CIGAR.append(str(g_value)+g_type)
				g_value=0
			elif g_type=='':
				output_CIGAR.append(str(t_value)+t_type)
				t_value=0
			else:
				print "BUG:"+transcript_CIGAR
				break
			
			if g_value<=0:
				pop_g_structure=True
			if t_value<=0:
				pop_t_structure=True
		
	except TypeError:
		print "TypeError:"
	
	return "".join(output_CIGAR)
	

def build_CIGAR(pos_in_exon,exon_number,DNA_structure,read_length=75,align_sense="+"):
	'''
	Compute the CIGAR string based on DNA_structure (a.k.a ENST_structures[ENST])
	Cf. http://samtools.sourceforge.net/SAM1.pdf
	
	Now we know this read is start from the i-th exon, and we know the read starting at the 
	position (1-based) in the i-th exon. 
	
	We always presume that the pos_in_exon is the position of the first base of the read when the gene is from postive strand, and the pos_in_exon is the position of the last base of the read when the gene is from reversed strand.
	'''
	global CIGAR_parser,Number_parser
	
	CIGAR=[]
	tmp=0
	tmp_bp=read_length
	exon_count=len(DNA_structure[1])
	if DNA_structure[2]=="+":
		# Compute the length of mathced subsequence in this exon.
		if exon_number==0:
			tmp=DNA_structure[1][0]-pos_in_exon+1
		else:
			tmp=(DNA_structure[1][exon_number]-DNA_structure[1][exon_number-1])-pos_in_exon+1
		
		if tmp<=0:
			# pos_in_exon is larger than the length of exon `exon_number`.
			return False
		
		for j in xrange(exon_number,exon_count):
			if tmp_bp<=tmp:
				CIGAR.append(str(tmp_bp)+"M")
				break
			else:
				CIGAR.append(str(tmp)+"M")
				tmp_bp=tmp_bp-tmp
				if j+1>=exon_count:
					return False
				else:
					# compute the length of the j-th intron.
					if j==0:
						tmp=DNA_structure[0][1]-(DNA_structure[0][0]+DNA_structure[1][0])
					else:
						tmp=DNA_structure[0][j+1]-(DNA_structure[0][j]+(DNA_structure[1][j]-DNA_structure[1][j-1]))
				CIGAR.append(str(tmp)+"N")
				
				# Compute the length of the next exon.
				tmp=(DNA_structure[1][j+1]-DNA_structure[1][j])
			
	else:
		'''
		Note that, since we trace current exon down to 0-th exon according postive strand position, the returned CIGAR format will coordinate to positive strand.
		'''
		# Compute the length of mathced subsequence in this exon.
		if exon_number==0:
			tmp=DNA_structure[1][0]-pos_in_exon+1
		else:
			tmp=(DNA_structure[1][exon_number]-DNA_structure[1][exon_number-1])-pos_in_exon+1
		
		if tmp<=0:
			# pos_in_exon is larger than the length of exon `exon_number`.
			return False
		
		for j in xrange(exon_number,-1,-1):
			#print tmp
			if tmp_bp<=tmp:
				CIGAR.append(str(tmp_bp)+"M")
				break
			else:
				CIGAR.append(str(tmp)+"M")
				tmp_bp=tmp_bp-tmp
				if j-1<0:
					return False
				else:
					# compute the length of the intron between j-th exon and (j-1)-th exon
					if j==0:
						# There is no intron between 0-th exon and (-1)-th exon
						return False
					else:
						tmp=DNA_structure[0][j-1]-(DNA_structure[0][j]+(DNA_structure[1][j]-DNA_structure[1][j-1]))
				CIGAR.append(str(tmp)+"N")
				
				# Compute the length of the next exon.
				if j-2>=0:
					tmp=(DNA_structure[1][j-1]-DNA_structure[1][j-2])
				else:
					tmp=DNA_structure[1][j-1]
	
	return "".join(CIGAR)
	


def genome_GPS(ENST,pos,ENST_structures,read_length=75,align_sense="+"):
	'''
	Given ENST name and position in ENST, this function will return the genomic position
	Note that pos is 1-based coordinates.
	'''
	if ENST not in ENST_structures:
		print ENST," not in Your annoation."
		#tmp=ENST_structures.keys()
		#print tmp[:100]
		return False
	else:
		'''
			ENST_structures: prev_ENST -> (exon_starts,exon_widths,prev_ENST_strand,prev_ENST_chromosme)
			Note that exon_widths is accumulated.
		'''
		if ENST_structures[ENST][2]=="-":
			pos=pos+read_length-1
		
		offset=0
		tmp_len=len(ENST_structures[ENST][1])
		for i in xrange(tmp_len):
			if pos<=ENST_structures[ENST][1][i]:
				'''
					Given a read ACTTCCCTTG, we say the rightest base is the first base of this read,
					and the leftest base is the last base of this read, i.e., A and G respectively.
					
					The SAM format will return the position of the first base in cDNA cooridnate.
					However, when we align this read to the cDNA sequence of a gene on reversed strand,
					we need to convert the position to the last base of this read to make the notations
					consistent to positive strand.
				'''
				if ENST_structures[ENST][2]=="+":
					if i==0:
						offset=pos-1
					else:
						offset=pos-ENST_structures[ENST][1][i-1]-1
				else:
					if i==0:
						offset=ENST_structures[ENST][1][0]-pos
					else:
						offset=ENST_structures[ENST][1][i]-pos
				pos_in_exon=offset+1
				return [ENST_structures[ENST][3],ENST_structures[ENST][0][i]+(offset),build_CIGAR(pos_in_exon,i,ENST_structures[ENST],read_length)]
		
	# This line will be executed only if pos > the length of whole tanscript ENST.
	return False


def verification_for_bowtie_PIPE(input_file,output_file,ENST_structures):
	global options
	sucessful_count=0
	prev_alignment=None
	curr_block=[]

	fhr=subprocess.Popen(input_file, shell=True,stdout=subprocess.PIPE).stdout
	fhw=open(output_file,"w+")
	for line in fhr:
		if Is_SAM_header(line):
			continue
		line_tmp=line.strip(" \r\n\t")
		line_tmp=line_tmp.split("\t")
		if line_tmp[2]=='*':
			continue
		if int(line_tmp[1]) & 4 ==4:
			continue
		
		'''
		Since the output of bowtie will be one read by one reads, 
		the alignments of one read will be put together as a cluster or a block.
		Thus we compute one block by one block.
		'''
		if prev_alignment!=line_tmp[0]:
			# Enter a new block
			# Dump previously block
			tmp_len=len(curr_block)
			pos_set=set()
			if tmp_len>0:
				read_name=set()
				for alignment in curr_block:
					read_name.add(alignment[0])
					tmp_read_length=adjust_length_by_CIGAR(alignment[5],len(alignment[9]))
					target=genome_GPS(alignment[2],int(alignment[3]),ENST_structures,tmp_read_length)
					if target==False or target[1]==False or target[2]==False:
						print "FAIL in genome_GPS:",alignment,target
						read_name=[]
						pos_set=[]
						break
					pos_set.add((target[0],target[1]))
				
				if len(read_name)>1:
					#If this program wrongly put the alignments of different reads into the same block, then we just exit!
					print "Alignments of each read should be put together.\n",read_name,curr_block,pos_set
					raise SystemExit
				
				if len(pos_set)==1:
					# All alignments has the same genomic position. 
					sucessful_count+=1
					alignment=curr_block[0]
					if ENST_structures[alignment[2]][2]=="-":
						# Since we align cDNA file, we need to convert it to either positive or negitive strand accodring to ENST.
						alignment[5]=reverse_CIGAR(alignment[5])
						if (int(alignment[1]) & 16)==16:
							alignment[1]= str(int(alignment[1]) - 16)
						else:
							alignment[1]= str(int(alignment[1]) | 16)
						alignment[9]=complement_sequence(alignment[9])
						alignment[10]=alignment[10][::-1]	# reverse the Phread Qualities

					alignment[5]=merge_CIGARs(alignment[5],target[2])
					alignment[2]=ENST_structures[alignment[2]][3]
					tmp=pos_set.pop()
					alignment[3]=str(tmp[1])
					fhw.write("\t".join(alignment)+"\n")
				else:
					print read_name,"::",str(len(curr_block)),"::",pos_set
			
			# initialize for new block
			prev_alignment=line_tmp[0]
			curr_block=[line_tmp]
		else:
			curr_block.append(line_tmp)
	#fhr.close()
	# Dump the last block
	tmp_len=len(curr_block)
	pos_set=set()
	if tmp_len>0:
		read_name=set()
		for alignment in curr_block:
			read_name.add(alignment[0])
			tmp_read_length=adjust_length_by_CIGAR(alignment[5],len(alignment[9]))
			target=genome_GPS(alignment[2],int(alignment[3]),ENST_structures,tmp_read_length)
			if target==False or target[1]==False or target[2]==False:
				print alignment,target
				continue
			pos_set.add((target[0],target[1]))
		
		if len(read_name)>1:
			#If this program wrongly put the alignments of different reads into the same block, then we just exit!
			print "Alignments of each read should be put together.\n",read_name,curr_block,pos_set
			raise SystemExit
		
		if len(pos_set)==1:
			# All alignments has the same genomic position. 
			sucessful_count+=1
			alignment=curr_block[0]
			if ENST_structures[alignment[2]][2]=="-":
				# Since we align cDNA file, we need to convert it to either positive or negitive strand accodring to ENST.
				alignment[5]=reverse_CIGAR(alignment[5])
				if (int(alignment[1]) & 16)==16:
					alignment[1]= str(int(alignment[1]) - 16)
				else:
					alignment[1]= str(int(alignment[1]) | 16)
				alignment[9]=complement_sequence(alignment[9])
				alignment[10]=alignment[10][::-1]	# reverse the Phread Qualities

			alignment[5]=merge_CIGARs(alignment[5],target[2])
			alignment[2]=ENST_structures[alignment[2]][3]
			tmp=pos_set.pop()
			alignment[3]=str(tmp[1])
			fhw.write("\t".join(alignment)+"\n")
		else:
			print read_name,"::",str(len(curr_block)),"::",pos_set
	fhw.close()
	
	print sucessful_count
	return True



def verification_for_bowtie(input_file,output_u_file,ENST_structures,output_m_file=stdout):
	sucessful_count=0
	prev_alignment=None
	curr_block=[]
	fhr=open(input_file,"r")
	fhw=open(output_u_file,"w+")
	if output_m_file!=stdout:
		fhw2=open(output_m_file,"w+")
	else:
		fhw2=stdout
	
	for line in fhr:
		if Is_SAM_header(line):
			continue
		line_tmp=line.strip(" \r\n\t")
		line_tmp=line_tmp.split("\t")
		if line_tmp[2]=='*':
			continue
		if int(line_tmp[1]) & 4 ==4:
			continue
		
		'''
		Since the output of bowtie will be one read by one reads, 
		the alignments of one read will be put together as a cluster or a block.
		Thus we compute one block by one block.
		'''
		if prev_alignment!=line_tmp[0]:
			# Enter a new block
			# Dump previously block
			tmp_len=len(curr_block)
			pos_set=set()
			if tmp_len>0:
				read_name=set()
				for alignment in curr_block:
					read_name.add(alignment[0])
					tmp_read_length=adjust_length_by_CIGAR(alignment[5],len(alignment[9]))
					target=genome_GPS(alignment[2],int(alignment[3]),ENST_structures,tmp_read_length)
					if target==False or target[1]==False or target[2]==False:
						print "FAIL in genome_GPS:",alignment,target
						read_name=[]
						pos_set=[]
						break
					pos_set.add((target[0],target[1]))
				
				if len(read_name)>1:
					#If this program wrongly put the alignments of different reads into the same block, then we just exit!
					print "Alignments of each read should be put together.\n",read_name,curr_block,pos_set
					raise SystemExit
				
				if len(pos_set)==1:
					# All alignments has the same genomic position. 
					sucessful_count+=1
					alignment=curr_block[0]
					if ENST_structures[alignment[2]][2]=="-":
						# Since we align cDNA file, we need to convert it to either positive or negitive strand accodring to ENST.
						alignment[5]=reverse_CIGAR(alignment[5])
						if (int(alignment[1]) & 16)==16:
							alignment[1]= str(int(alignment[1]) - 16)
						else:
							alignment[1]= str(int(alignment[1]) | 16)
						alignment[9]=complement_sequence(alignment[9])
						alignment[10]=alignment[10][::-1]	# reverse the Phread Qualities

					alignment[5]=merge_CIGARs(alignment[5],target[2])
					isofrom_strand=ENST_structures[alignment[2]][2]
					alignment[2]=ENST_structures[alignment[2]][3]
					tmp=pos_set.pop()
					alignment[3]=str(tmp[1])
					if -1!=alignment[5].find("N"):
						alignment.append("XS:A:"+isofrom_strand)
					fhw.write("\t".join(alignment)+"\n")
				else:
					print >>fhw2, read_name,"::",str(len(curr_block)),"::",pos_set
			
			# initialize for new block
			prev_alignment=line_tmp[0]
			curr_block=[line_tmp]
		else:
			curr_block.append(line_tmp)
	fhr.close()
	# Dump the last block
	tmp_len=len(curr_block)
	pos_set=set()
	if tmp_len>0:
		read_name=set()
		for alignment in curr_block:
			read_name.add(alignment[0])
			tmp_read_length=adjust_length_by_CIGAR(alignment[5],len(alignment[9]))
			target=genome_GPS(alignment[2],int(alignment[3]),ENST_structures,tmp_read_length)
			if target==False or target[1]==False or target[2]==False:
				print alignment,target
				continue
			pos_set.add((target[0],target[1]))
		
		if len(read_name)>1:
			#If this program wrongly put the alignments of different reads into the same block, then we just exit!
			print "Alignments of each read should be put together.\n",read_name,curr_block,pos_set
			raise SystemExit
		
		if len(pos_set)==1:
			# All alignments has the same genomic position. 
			sucessful_count+=1
			alignment=curr_block[0]
			if ENST_structures[alignment[2]][2]=="-":
				# Since we align cDNA file, we need to convert it to either positive or negitive strand accodring to ENST.
				alignment[5]=reverse_CIGAR(alignment[5])
				if (int(alignment[1]) & 16)==16:
					alignment[1]= str(int(alignment[1]) - 16)
				else:
					alignment[1]= str(int(alignment[1]) | 16)
				alignment[9]=complement_sequence(alignment[9])
				alignment[10]=alignment[10][::-1]	# reverse the Phread Qualities

			alignment[5]=merge_CIGARs(alignment[5],target[2])
			isofrom_strand=ENST_structures[alignment[2]][2]
			alignment[2]=ENST_structures[alignment[2]][3]
			tmp=pos_set.pop()
			alignment[3]=str(tmp[1])
			if -1!=alignment[5].find("N"):
				alignment.append("XS:A:"+isofrom_strand)
			fhw.write("\t".join(alignment)+"\n")
		else:
			print >>fhw2, read_name,"::",str(len(curr_block)),"::",pos_set
	fhw.close()
	if fhw2!=stdout:
		fhw2.close()
	print sucessful_count
	return True


def adjust_alignment_for_strands(ENST_structures,alignment,target):
	if ENST_structures[alignment[2]][2]=="-":
		# Since we align cDNA file, we need to convert it to either positive or negitive strand accodring to ENST.
		alignment[5]=reverse_CIGAR(alignment[5])
		if (int(alignment[1]) & 16)==16:
			alignment[1]= str(int(alignment[1]) - 16)
		else:
			alignment[1]= str(int(alignment[1]) | 16)
		alignment[9]=complement_sequence(alignment[9])
		alignment[10]=alignment[10][::-1]	# reverse the Phread Qualities

	alignment[5]=merge_CIGARs(alignment[5],target[2])
	isofrom_strand=ENST_structures[alignment[2]][2]
	alignment[2]=ENST_structures[alignment[2]][3]
	alignment[3]=str(target[1])
	if -1!=alignment[5].find("N"):
		alignment.append("XS:A:"+isofrom_strand)
	
	return "\t".join(alignment)


def adjust_length_by_CIGAR(CIGAR,read_length):
	global Insertion_parser,Deletion_parser,Skip_parser,SoftClip_parser
	# Insertion will shorten the real mapped length of reads
	ip=Insertion_parser.findall(CIGAR)
	tmp_read_length=read_length
	if len(ip)!=0:
		for i in ip:
			tmp_read_length-=int(i)
	
	# Soft Clip will extend the real mapped length of reads
	ip=SoftClip_parser.findall(CIGAR)
	if len(ip)!=0:
		for i in ip:
			tmp_read_length-=int(i)
	
	# Deletion will extend the real mapped length of reads
	ip=Deletion_parser.findall(CIGAR)
	if len(ip)!=0:
		for i in ip:
			tmp_read_length+=int(i)
	
	# Skip will extend the real mapped length of reads
	ip=Skip_parser.findall(CIGAR)
	if len(ip)!=0:
		for i in ip:
			tmp_read_length+=int(i)
	
	return tmp_read_length
	

def verification_for_BWA(input_file,output_file,ENST_structures,convert_unique=False):
	global options,CIGAR_parser,SAM_XA_parser
	sucessful_count=0
	pattern=SAM_XA_parser
	
	fhr=open(input_file,"r")
	fhw=open(output_file,"w+")
	for line in fhr:
		if Is_SAM_header(line):
			continue
		
		line_tmp=line.strip(" \r\n\t")
		line_tmp=line_tmp.split("\t")
		
		if line_tmp[2]=='*':
			continue
		
		flag=int(line_tmp[1])
		if (flag>>2)%2==1:
			continue
		
		MRNM=line_tmp[2]
		MRNM_pos=int(line_tmp[3])
		read_length=len(line_tmp[9])
		
		if (-1==line.find("XA:Z")):
			'''
			There are two cases, 
			The first one is there are too much targets. We will skip this line.
			The second one is this one is unique alignmentm. Sincere we assume 
			the input file are only reads with multiple targets, we skip this line, too.
			'''
			continue
		else:
			#m=pattern.findall(line)
			m=extract_SAM_XA(line)
			if len(m)==0:
				print "WRONG FORMAT:",line
				continue
			else:
				pos_set=set()
				
				tmp_read_length=adjust_length_by_CIGAR(line_tmp[5],read_length)
				major_target=genome_GPS(MRNM,abs(int(MRNM_pos)),ENST_structures,tmp_read_length)
				if major_target==False or major_target[1]==False or major_target[2]==False:
					print line_tmp,major_target
					continue
				pos_set.add(tuple(major_target[:2]))
				
				for i in m:
					i=i.split(",")
					tmp_read_length=adjust_length_by_CIGAR(i[2],read_length)
					tmp=genome_GPS(i[0],abs(int(i[1])),ENST_structures,tmp_read_length)
					if tmp==False or tmp[1]==False or tmp[2]==False:
						break
					pos_set.add(tuple(tmp[:2]))
				
				# Generate Output
				if len(pos_set)==1:
					sucessful_count=sucessful_count+1
					if ENST_structures[MRNM][2]=="-":
						# Since we align cDNA file, we need to convert it to either positive or negitive strand accodring to ENST.
						line_tmp[5]=CIGAR_parser.findall(line_tmp[5])
						line_tmp[5].reverse()
						line_tmp[5]="".join(line_tmp[5])
						if (int(line_tmp[1]) & 16)==16:
							line_tmp[1]= str(int(line_tmp[1]) - 16)
						else:
							line_tmp[1]= str(int(line_tmp[1]) | 16)
						
						line_tmp[9]=list(line_tmp[9])
						line_tmp[9].reverse()
						for i in xrange(len(line_tmp[9])):
							if line_tmp[9][i]=="N":
								continue
							elif line_tmp[9][i]=="A":
								line_tmp[9][i]="T"
							elif line_tmp[9][i]=="G":
								line_tmp[9][i]="C"
							elif line_tmp[9][i]=="C":
								line_tmp[9][i]="G"
							elif line_tmp[9][i]=="T":
								line_tmp[9][i]="A"
							else:
								print "Wrong Sequence Format"
								exit(2)
						line_tmp[9]="".join(line_tmp[9])
						line_tmp[10]=line_tmp[10][::-1]	# reverse the Phread Qualities
					
					line_tmp[5]=merge_CIGARs(line_tmp[5],major_target[2])
					line_tmp[2]=ENST_structures[line_tmp[2]][3]
					line_tmp[3]=str(major_target[1])
					
					fhw.write("\t".join(line_tmp)+"\n")
				else:
					if options.quite==False:
						print MRNM,"::",len(m),"::",pos_set
					
				
	fhw.close()
	fhr.close()
	
	print sucessful_count
	return True


def build_G2T_structures(ENST_structures):
	intervals=dict()
	for ENST_id in ENST_structures.iterkeys():
		ENST=ENST_structures[ENST_id]
		if ENST[3] not in intervals:
			intervals[ENST[3]]=[]

		for i in xrange(len(ENST[0])):
			if i==0:
				intervals[ENST[3]].append(Interval(ENST[0][0],ENST[0][0]+ENST[1][0]-1,ENST_id+":0"))
			else:
				intervals[ENST[3]].append(Interval(ENST[0][i],ENST[0][i]+(ENST[1][i]-ENST[1][i-1])-1,ENST_id+":"+str(i)))	
	
	G2T_tree=dict()
	for i in intervals.iterkeys():
		G2T_tree[i]=IntervalTree(intervals[i],24)
	
	del intervals
	return G2T_tree


def exon_over_genome_position(pos,G2T_tree,ENST_structures,offset_to_pos,chr_name=""):
	output=[]
	if chr_name=="":
		result=[]
		for i in G2T_tree.iterkeys():
			result+=G2T_tree[i].find(pos,pos+offset_to_pos)
	else:
		if chr_name not in G2T_tree:
			return output
		result=G2T_tree[chr_name].find(pos,pos+offset_to_pos)
	for i in result:
		output.append((i.id).split(":"))
	return output


def convert_genome_position_to_transcirpt_position(pos,chr_name,G2T_tree,ENST_structures,read_end_offset=0):
	'''
	read_end_offset == (read_length-1)
	'''
	output=[]
	result=G2T_tree[chr_name].find(pos,pos+read_end_offset)
	for i in result:
		tmp=(i.id).split(":")
		j=int(tmp[1])
		ENST=ENST_structures[tmp[0]]
		
		if ENST[2]=="+":
			output.append((tmp[0],ENST[1][j]-(i.stop-pos)))
		if ENST[2]=="-":
			output.append((tmp[0],ENST[1][j]-(pos+read_end_offset-i.start)))
	return output


def test_G2T(ENST_structures):
	from random import randint,choice
	G2T_tree=build_G2T_structures(ENST_structures)
	chr_set=G2T_tree.keys()
	for i in xrange(10):
		chr_name=choice(chr_set)
		test_pos=[randint(1,1000000) for i in xrange(1000)]
		for pos in test_pos:
			tmp=convert_genome_position_to_transcirpt_position(pos,chr_name,G2T_tree,ENST_structures)
			try:
				for j in tmp:
					tmp2=genome_GPS(j[0],j[1],ENST_structures)
					if tmp2==False:
						print "X",pos,j,ENST_structures[j[0]][2]
						continue
					g,c=tmp2
					if c==False:
						print "X",pos,j,ENST_structures[j[0]][2],g,c
						continue
					if g!=pos:
						return False
						print "X",pos,j,ENST_structures[j[0]][2],g
					else:
						print pos,j,ENST_structures[j[0]][2],g
			except TypeError:
				print "XX",j 
	return True


if __name__ == '__main__':
	usage = "%prog [options] input.sam gene_annoation.gtf output.sam"
	parser = OptionParser(usage=usage)
	parser.add_option("-s", "--source_type", dest="source_type", help="Specify the source of input file, 'BWA' or 'Bowtie'.", metavar="TYPE",default="Bowtie")
	parser.add_option("-b", "--annotation_type", dest="annotation_type", help="Specify the format of annotation file, 'BED' or 'GTF'.", metavar="TYPE",default="GTF")
	parser.add_option("-t", "--convert_unique", dest="convert_unique", action="store_true",help="Whether to convert positions of unique reads to genomic position. We assume alignment without 'XA:Z' as unique alignments, however it could be the case that a read has too much alignments.",default=False)
	parser.add_option("-q", "--quite", dest="quite", action="store_true",help="Make the program not print non-unique algined reads to stdout.",default=False)
	(options, args) = parser.parse_args()

	output_file=None
	filelist=[]
	if len(args) < 3:
		parser.error("missing required arguments.")
		exit(2)
	elif len(args)==3:
		input_file = args[0]
		gtf_file = args[1]
		output_file = args[2]
	else:
		parser.error("Too much arguments.")
		exit(2)
	
	if options.annotation_type=="GTF":
		ENST_structures=extract_gene_structure(gtf_file)
	elif options.annotation_type=="BED":
		ENST_structures=extract_gene_structure(gtf_file,"BED")
	else:
		parser.error("Wrong annotation format.")
		exit(2)
	print "Finished loading gene structure"
	
	print "Source Type:",options.source_type
	if options.source_type=="BWA":
		verification_for_BWA(input_file,output_file,ENST_structures,options.convert_unique)
	else:
		verification_for_bowtie(input_file,output_file,ENST_structures)
	
	
