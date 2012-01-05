'''
This program is for convert BLAT's result into standard SAM.
Thus it must be run after finished filtering



'''
import re
from optparse import OptionParser
from selected_the_multiply_mapped_reads_of_the_same_genomic_position import complement_sequence
CIGAR_parser = re.compile(r'([0-9]+[M|I|D|N|S|H|P])')
H_parser = re.compile(r'([0-9]+)H')
head_H = re.compile(r'^([0-9]+)H')


if __name__ == '__main__':
	usage = "%prog [options] input.sam source.fastq"
	parser = OptionParser(usage=usage)
	parser.add_option("-s", "--source_type", dest="source_type", help="Specify the source of input file, 'BWA' or 'Bowtie'.", metavar="TYPE",default="Bowtie")
	(options, args) = parser.parse_args()
	if len(args) < 2:
		parser.error("missing required arguments.")
		exit(2)
	elif len(args)==2:
		input_file = args[0]
		source_fq= args[1]

	if options.source_type!="Bowtie":
		print "Not yet implemented"
		exit(2)

	f=open(input_file,"r")
	reads_id=set()
	for line in f:
		reads_id.add(line[:line.find("\t")])
	f.close()

	f=open(source_fq,"r")
	FETCH=None
	reads_info=dict()
	l_count=0
	for line in f:
		l_count+=1
		if (l_count%4==1):
			id=line[1:].strip("\r\n")
			if id in reads_id:
				reads_info[id]=[]
				FETCH=id
		elif (l_count%4==2) and FETCH!=None:
			reads_info[FETCH].append(line.strip("\r\n"))
		elif (l_count%4==0) and FETCH!=None:
			reads_info[FETCH].append(line.strip("\r\n"))
			FETCH=None
	f.close()
	del reads_id

	f=open(input_file,"r")
	for line in f:
		line=line.strip(" \r\n")
		li=line.split("\t")
		li[4]="255"			# MAPQ field will affect the SNP calinng. Assign 255 to it for disabling this function.
		if int(li[1]) & 16 == 16:
			li[9]=complement_sequence(reads_info[li[0]][0])
			li[10]=reads_info[li[0]][1][::-1]
		else:
			li[9],li[10]=reads_info[li[0]]
		
		Hs=H_parser.findall(li[5])
		if Hs==[]:
			li[5]=li[5].replace('D','N')
			print "\t".join(li)
			continue
		count=0
		for H in Hs:
			count+=int(H)

		for i in xrange(len(li[11:])):
			i+=11
			tmp=li[i].find("NM:i:")
			if tmp!=-1:
				li[i]="NM:i:"+str(int(li[i][tmp+5:])+count)
				break
		
		tmp=head_H.findall(li[5])
		if tmp!=[]:
			li[3]=str(int(li[3])-int(tmp[0]))
			if li[3]<=0:
				continue

		li[5]=CIGAR_parser.findall(li[5].replace('H','M').replace('D','N'))
		for i in xrange(1,len(li[5])):
			first,second=li[5][i-1:i+1]
			if first[-1]==second[-1] and first[-1]=="M":
				li[5][i]=str(int(first[:-1])+int(second[:-1]))+'M'
				li[5][i-1]=""
		li[5]="".join(li[5])
		print "\t".join(li)
	f.close()

