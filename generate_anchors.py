'''
This program is to generate anchors for stage 2 of our algorithm.
This program splices strings into multiple substrings according to the given anchor length
'''
import subprocess
from optparse import OptionParser

if __name__ == '__main__':
	usage = "%prog [options] input_file output_file"
	parser = OptionParser(usage=usage)
	parser.add_option("-i", "--input_file_type", dest="input_file_type", help="FASTA/FASTQ/SAM [SAM]", metavar="TYPE",default="SAM")
	parser.add_option("-o", "--output_file_type", dest="output_file_type", help="FASTA/FASTQ [FASTQ]", metavar="TYPE",default="FASTQ")
	parser.add_option("-a", "--anchor_length", dest="anchor_length", help="anchor_length [25]", metavar="INT",default=25)
	parser.add_option("-r", "--read_length", dest="read_length", help="read_length [75]", metavar="INT",default=75)
	parser.add_option("-m", "--mismatch", dest="mismatch", help="mismatch [3]", metavar="INT",default=3)
	(options, args) = parser.parse_args()

	output_file=None
	snp_files=[]
	if len(args) != 2:
		parser.error("missing required arguments.")
		exit(2)
	elif len(args)==2:
		input_file = args[0] 
		output_file = args[1] 
	read_length=int(options.read_length)
	anchor_length=int(options.anchor_length)
	mismatch=int(options.mismatch)

	anchors=[]
	tmp=0
	while (tmp<=read_length):
		if tmp+anchor_length>read_length:
			break
		anchors.append((tmp,tmp+anchor_length))
		tmp+=anchor_length
	
	if options.input_file_type=="SAM":
		fw=open(output_file,"w+")
		instrcution="""awk '/^[^@]/{if (rshift($2, 2)%2==1){print $1;print $10; print $11; }}' """+input_file
		ps1 = subprocess.Popen(instrcution, shell=True,stdout=subprocess.PIPE).stdout
		count=0
		lines=[]
		for line in ps1:
			count+=1
			lines.append(line.strip(" \r\n"))
			if count %3==0:
				tmp=0
				head="@"+lines[0]+"__"
				for i,j in anchors:
					tmp+=1
					fw.write(head+str(tmp)+"\n"+lines[1][i:j]+"\n+\n"+lines[2][i:j]+"\n")
				lines=[]
				count=0
		fw.close()
	elif options.input_file_type=="FASTQ":
		fw=open(output_file,"w+")
		fr=open(input_file,"r")
		block=[]
		for line in fr:
			block.append(line.strip("\r\n"))
			if len(block)==4:
				if block[1].count("N")>mismatch:
					block=[]
					continue
				head=block[0].strip("\r\n") + "__"
				tmp=0
				for i,j in anchors:
					tmp+=1
					fw.write(head+str(tmp)+"\n"+block[1][i:j]+"\n+\n"+block[3][i:j]+"\n")
				block=[]
		fr.close()
		fw.close()




