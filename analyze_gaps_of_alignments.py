from sys import stdout , stderr
from math import log
from re import *
from optparse import OptionParser

pattern_A=compile(r"([0-9]+)M")
pattern_B=compile(r"([0-9]+)I")
pattern_C=compile(r"([0-9]+)[D|N]")

def CIGAR_operation_distribution_report(distribution):
	matches=distribution.keys()
	matches.sort()
	
	for i in matches:
		print i,"matches"
		
		intervals=distribution[i].keys()
		intervals.sort()
		
		for j in intervals:
			print int(pow(2,j)),",",distribution[i][j]
		
		print "total",sum(distribution[i].values())
		
		print "\nAccumlated"
		
		acc=0
		intervals.reverse()
		for j in intervals:
			acc+=distribution[i][j]
			print int(pow(2,j)),",",acc
		
		print "\n"
	

if __name__ == '__main__':
	usage = "%prog [options] input_file.sam \n"
	parser = OptionParser(usage=usage)
	
	#parser.add_option("-s", "--source_type", dest="source_type", help="Specify the source of input file, 'Bowtie' or 'BWA'.", metavar="TYPE",default="BWA")
	
	parser.add_option("-m", "--threshold_of_matches", dest="threshold_of_matches", help="Specify the threshold of matches. We only consider the reads which satisified this criteria. (>)", metavar="INT",default=73)
	parser.add_option("-S", "--run_statistics",action="store_true", dest="run_statistics", help="run and print the statistics of alignment", default=False)
	parser.add_option("-f", "--filter_unqualified_alignments", action="store_true",dest="filter_unqualified_alignments", help="filter unqualified alignments", default=False)
	parser.add_option("-g", "--min_gap_length", dest="min_gap_length", help="This parameter works only when -f is set true. If a read with a gap not less than min_gap_length bps, then we will print it (>=)[default=10]", default=10)
	(options, args) = parser.parse_args()
	
	output_file=None
	if len(args) < 1:
		parser.error("missing required arguments.")
		exit(2)
	elif len(args)==1:
		input_file = args[0]
	else:
		parser.error("Too much arguments.")
		exit(2)
	
	if (options.run_statistics==True) and (options.filter_unqualified_alignments==True):
		print >> stderr, "Parameters -S and -f are mutually excluded"
		exit(2)
	elif (options.run_statistics==True):
		print >> stderr, "Running statistics"
	elif (options.filter_unqualified_alignments==True):
		print >> stderr, "Filtering unqualified alignments"
	elif (options.run_statistics==False) and (options.filter_unqualified_alignments==False):
		print >> stderr, "At least one of parameters -S and -f should be set true."
		exit(2)
	
	print >> stderr, "options.threshold_of_matches:",options.threshold_of_matches
	print >> stderr, "options.min_gap_length:",options.min_gap_length
	threshold_of_matches=int(options.threshold_of_matches)
	min_gap_length=int(options.min_gap_length)
	
	distribution_B=dict()
	distribution_C=dict()
	
	f=open(input_file,"r")
	for li in f:
		li=li.strip(" \r\n")
		l=li.split("\t")
		l=l[5]
		m=pattern_A.findall(l)
		if len(m)==0:
			continue
		
		m_c=0
		for i in m:
			m_c+=int(i)
		
		if options.filter_unqualified_alignments==True:
			if m_c > threshold_of_matches:
				m=pattern_C.findall(l)
				gap=[0]
				for i in m:
					gap.append(int(i))
				if max(gap)>=min_gap_length:
					print li
			continue
			
		
		if m_c > threshold_of_matches:
			if m_c not in distribution_B:
				distribution_B[m_c]=dict()
			if m_c not in distribution_C:
				distribution_C[m_c]=dict()
			
			
			# Compute Pattern_B
			m=pattern_B.findall(l)
			if len(m)!=0:
				gap=[]
				for i in m:
					gap.append(int(i))
				
				gap=int(log(max(gap),2))
				if gap not in distribution_B[m_c]:
					distribution_B[m_c][gap]=0
				distribution_B[m_c][gap]+=1
			
			# Compute Pattern_C
			m=pattern_C.findall(l)
			if len(m)!=0:
				gap=[]
				for i in m:
					gap.append(int(i))
				gap=int(log(max(gap),2))
				if gap not in distribution_C[m_c]:
					distribution_C[m_c][gap]=0
				distribution_C[m_c][gap]+=1
		
	f.close()
	if (options.run_statistics==True):
		print "Insertion Distribution"
		CIGAR_operation_distribution_report(distribution_B)
		print "Deletion/Skip Distribution"
		CIGAR_operation_distribution_report(distribution_C)

