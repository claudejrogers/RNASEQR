import commands
import subprocess
import re
from optparse import OptionParser


def filter_and_write(input_file,output_file,reads_count,condition_func):
	try:
		f=open(input_file,'r')
		fw=open(output_file,'w+')
	except TypeError:
		print input_file
	c=0
	for line in f:
		line_tmp=line.strip(" \r\n\t")
		line_tmp=line_tmp.split("\t")
		
		if line_tmp[0] in reads_count:
			if condition_func(reads_count[line_tmp[0]]):
				c+=1
				fw.write(line)
	
	print "Counts:",c
	f.close()
	fw.close()
	return True
	
def report_multi_reads(input_file,reads_count):
	tmp=input_file.split(".")
	tmp.append("multiple."+tmp.pop())
	output_file=".".join(tmp)
	return filter_and_write(input_file,output_file,reads_count,lambda x: x>1)

def report_unique_lines(input_file,reads_count):
	tmp=input_file.split(".")
	tmp.append("unique."+tmp.pop())
	output_file=".".join(tmp)
	return filter_and_write(input_file,output_file,reads_count,lambda x: x==1)

def count_unique_reads(reads_count):
	unique=0
	for i in reads_count.iterkeys():
		if reads_count[i]==1:
			unique=unique+1
	return unique

def Is_SAM_header(line,extension=[]):
	SAM_spec=["@HD","@SQ","@RG","@PG","@CO"]
	SAM_spec=SAM_spec+extension
	
	if line[:3] in SAM_spec:
		return True
	return False


if __name__ == '__main__':
	usage = "%prog [options] <input_prefix.sam> "
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	
	output_file=None
	filelist=[]
	if len(args) < 1:
		parser.error("missing required arguments.")
		exit(2)
	elif len(args)==1:
		input_file = args[0]
	else:
		parser.error("Too much arguments.")
		exit(2)
	
	if True:
		tmp=input_file.split(".")
		tmp.append("unique."+tmp.pop())
		output_file_u=".".join(tmp)
		fw_u=open(output_file_u,'w+')
		
		c=0
		f = open(input_file,"r")
		prev=""
		curr_block=[]
		for line in f:
			if (Is_SAM_header(line)):
				continue
			li=line.split("\t")
			if int(li[1]) & 4 ==4:
				continue
			head=li[0]
			if prev!=head:
				if len(curr_block)==1:
					fw_u.write(curr_block[0])
				prev=head
				curr_block=[line]
			else:
				curr_block.append(line)
		if len(curr_block)==1:
			fw_u.write(curr_block[0])
		fw_u.close()
	
