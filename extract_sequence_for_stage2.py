'''
This program is used to develope the methods of stage2.
It need to be generalized for further applications,e.g., more than three segments. 
'''
from optparse import OptionParser
import subprocess
from sys  import stderr,argv
from map_and_tables import *
from random import randint

def processing_bundle(bundle):
	if len(bundle)==1:
		for k in bundle.iterkeys():
			lines=bundle[k]
			bucket=set()
			for li in lines:
				bucket.add((li[1],li[2]))
			if len(bucket)==1:
				for i in xrange(len(lines)):
					lines[i]="\t".join(lines[i])
				return lines
	else:
		print "ERROR",bundle
	return False


def generate_fasta(data_file,reads,out_fasta):
	f=open(data_file,"r")
	fw=open(out_fasta,"w+")
	count=0
	bucket=[]
	for l in f:
		bucket.append(l)
		if count%4==3:
			x=bucket[0][1:].strip("\r\n ")
			if x in reads:
				fw.write(">"+bucket[0][1:])
				fw.write(bucket[1])
			bucket=[]
		count+=1
	f.close()
	fw.close()
	return True


def generate_fastas(data_file,reads,out_fastas):
	f=open(data_file,"r")
	count=0
	bucket=[]
	for l in f:
		bucket.append(l)
		if count%4==3:
			x=bucket[0][1:].strip("\r\n ")
			if x in reads:
				out_fastas[reads[x]].write(">"+bucket[0][1:]+bucket[1])
			bucket=[]
		count+=1
	f.close()
	return True


def calculate_distribution_of_difference(filename):
	dist=dict()
	s2g_reads=dict()
	s2g_reads_ori=dict()
	#instrcution="""cut -f 1,4 """+output_file[:output_file.rfind(".")]+".2g.sam"
	instrcution="""cut -f 1,2,4 """+filename
	ps3 = subprocess.Popen(instrcution, shell=True,stdout=subprocess.PIPE).stdout
	for line in ps3:
		f=line.strip(" \r\n\t").split("\t")
		x=f[0].split("__")
		name=x[0]
		seq_no=int(x[1])
		if name not in s2g_reads:
			s2g_reads[name]=dict()
			s2g_reads_ori[name]=set()
		s2g_reads[name][seq_no]=int(f[2])
		s2g_reads_ori[name].add(f[1])
	# Validate
	for i in s2g_reads_ori.iterkeys():
		if len(s2g_reads_ori[i])>1:
			#print >>stderr, "ERROR_ORI",s2g_reads_ori[i],s2g_reads[i]
			del s2g_reads[i]
		else:
			s2g_reads_ori[i]=int(s2g_reads_ori[i].pop())
	#print >>stderr, len(s2g_reads)
	
	#count
	counts=0
	f_counts=0
	for y in s2g_reads.iterkeys():
		x=s2g_reads[y]
		k=x.keys()
		k.sort()
		for i in xrange(len(k)-1,0,-1):
			diff=x[k[i]]-x[k[i-1]]
			if (diff<0 and s2g_reads_ori[y]!=16) or (diff>0 and s2g_reads_ori[y]!=0):
				#print >>stderr, s2g_reads[y],s2g_reads_ori[y]
				f_counts+=1
				continue
			if diff not in dist:
				dist[diff]=0
			dist[diff]+=1
			counts+=1
	del s2g_reads
	print "FAIL ORI.,",f_counts
	interval=counts/100
	index=sorted(dist.keys())
	counts=1
	for i in xrange(1,len(index)):
		dist[index[i]]+=dist[index[i-1]]
	print "Interval",interval
	print "quantile,value,acc"
	i=0
	while (i<len(index)):
		if (counts*interval)<=dist[index[i]]:
			print counts,",",index[i],",",counts*interval
			counts+=1
		else:
			i+=1
	return True


def from_T_or_from_G(phase3_input_file,test_file):
	phase3_result=dict()

	# Get the valid reafs id
	instrcution="""cut -f 1,2,3 """+phase3_input_file
	ps1 = subprocess.Popen(instrcution, shell=True,stdout=subprocess.PIPE).stdout
	for line in ps1:
		li=line.strip(" \r\n").split("\t")
		tmp=li[0].split('__')
		if tmp[0] not in phase3_result:
			phase3_result[tmp[0]]=[]
		phase3_result[tmp[0]].append((li[1],li[2]))
	tmp=phase3_result.keys()
	for i in tmp:
		len_l=len(phase3_result[i])
		len_s=len(set(phase3_result[i]))
		if len_s!=1:
			del phase3_result[i]
		else:
			phase3_result[i]=None
	
	# Get the valid segment IDs
	valid_segment_id_from_phase3_input=set()
	instrcution="""cut -f 1 """+phase3_input_file
	ps1 = subprocess.Popen(instrcution, shell=True,stdout=subprocess.PIPE).stdout
	for line in ps1:
		li=line.split("__")
		if li[0] in phase3_result:
			valid_segment_id_from_phase3_input.add(line.strip("\r\n\t "))
	del phase3_result
	
	dist=dict(zip(["TT","GG","GT","TG"],[0,0,0,0]))
	reads=dict()
	instrcution="""cut -f 1 """+test_file
	ps1 = subprocess.Popen(instrcution, shell=True,stdout=subprocess.PIPE).stdout
	for line in ps1:
		line=line.strip(" \r\n")
		li=line.split('__')
		li=li[0]
		if li not in reads:
			reads[li]=[]
		if line in valid_segment_id_from_phase3_input:
			reads[li].append("T")
		else:
			reads[li].append("G")
	for i in reads.iterkeys():
		x="".join(reads[i])
		dist[x]+=1
	print "TT,",dist["TT"]
	print "GG,",dist["GG"]
	print "GT/TG,",dist["GT"]+dist["TG"]
	return True	




if __name__ == '__main__':
	usage = "%prog [options] source.fastq anchors_on_transcriptome.sam anchors_on_genome.sam output_prefix "
	parser = OptionParser(usage=usage)
	parser.add_option("-p", "--path", dest="local_align_path", help="local_align_path", metavar="PATH",default="/local/blat_v34/blat")
	parser.add_option("-i", "--index", dest="local_align_index", help="local_align_index", metavar="PATH",default="/local/seq_indexes/blat_indexes/Homo_sapiens.GRCh37.59.genome.fa.2bit")
	parser.add_option("-n", "--num_seq", dest="num_of_seg", help="The number of valid segments. It will be used to filter anchors", metavar="INT",default=2)
	parser.add_option("-d", "--debug", dest="debug", help="", metavar="INT",default=-1)
	(options, args) = parser.parse_args()

	if len(args) < 4:
		parser.error("missing required arguments.")
		exit(2)
	elif len(args)==4:
		source_file=args[0]
		phase3_input_file = args[1]
		phase4_input_file = args[2]
		output_file = args[3]
		num_of_seg = int(options.num_of_seg)
		debug = int(options.debug)
	else:
		parser.error("Too much arguments.")
		exit(2)
	

	
	#
	# Fetch the IDs of the vaild reads from  anchors_on_transcriptome.sam
	#
	phase3_result=dict()
	instrcution="""cut -f 1,2,3 """+phase3_input_file
	ps1 = subprocess.Popen(instrcution, shell=True,stdout=subprocess.PIPE).stdout
	count_1,count_2,count_3=(0,0,0)
	f_count_1,f_count_2,f_count_3=(0,0,0)
	for line in ps1:
		li=line.strip(" \r\n").split("\t")
		tmp=li[0].split('__')
		if tmp[0] not in phase3_result:
			phase3_result[tmp[0]]=[[],[]]	# the first list is for storing orientation and ref. seq. ID. The second list is to store the segment ID.
		phase3_result[tmp[0]][0].append((li[1],li[2]))
		phase3_result[tmp[0]][1].append(tmp[1])
	# To filter the reads whosesegments are not mapped to the same reference with the same orination.
	for i in phase3_result.iterkeys():
		len_l=len(phase3_result[i][0])
		len_s=len(set(phase3_result[i][0]))
		if len_l==1:
			pass
		elif len_s==1:
			pass
		else:
			phase3_result[i]=[[],[]]

	#
	# Use the filtered read to fetch the candiantes segments to outputfile
	#
	candidate_reads=[]
	# scan phase4 results at first.
	bundle=dict()
	f=open(phase4_input_file,"r")
	for line in f:
		li=line.strip(" \r\n").split("\t")
		tmp=li[0].find('__')
		if (li[0][:tmp] in phase3_result):
			if (li[0][tmp+2:] in phase3_result[li[0][:tmp]][1]):
				continue
		if li[0][:tmp] not in bundle:
			tmp2=processing_bundle(bundle)
			if tmp2 != False:
				for tmp3 in tmp2:
					candidate_reads.append(tmp3)
			bundle.clear()
			bundle[li[0][:tmp]]=[]
		bundle[li[0][:tmp]].append(li)
	if  len(bundle)!=0:
		tmp2=processing_bundle(bundle)
		if tmp2 != False:
			for tmp3 in tmp2:
				candidate_reads.append(tmp3)
		bundle.clear()
	f.close()

	# scan pase3 result now
	bundle=dict()
	f=open(phase3_input_file,"r")
	for line in f:
		li=line.strip(" \r\n").split("\t")
		tmp=li[0].find('__')
		if li[0][:tmp] not in bundle:
			tmp2=processing_bundle(bundle)
			if tmp2 != False:
				for tmp3 in tmp2:
					candidate_reads.append(tmp3)
			bundle.clear()
			bundle[li[0][:tmp]]=[]
		bundle[li[0][:tmp]].append(li)
	if  len(bundle)!=0:
		tmp2=processing_bundle(bundle)
		if tmp2 != False:
			for tmp3 in tmp2:
				candidate_reads.append(tmp3)
		bundle.clear()
	f.close()
	del phase3_result


	# The selected candiantes segments are possbily containg the meaningless combinations.
	# If one of two segments of a read is from phase 3 and the other is from phase 4,
	# then, the combiantion could be mapped to dfferent reference sequences, or different 
	# orientatation.
	curr_result=dict()
	for li in candidate_reads:
		li=li.split("\t")
		lix=li[0][:li[0].find("__")]
		if lix not in curr_result:
			curr_result[lix]=[]
		curr_result[lix].append((li[0],li[1],li[2],li[3]))
	del candidate_reads
	tmp=curr_result.keys()
	print "CURRENT READS,",len(tmp)
	f_count_1=0
	for i in tmp:
		u=set()
		for x in curr_result[i]:
			u.add((x[1],x[2]))
		if len(u)!=1:
			del curr_result[i]
			f_count_1+=1
	del tmp
	print "DELETED READS,",f_count_1
	print "REMAIN READS,",len(curr_result)











	rand_sig=str(randint(1,1000000000))
	file_ptr=dict()
	for name in chr_name_id_map.itervalues():
		file_ptr[name]=open("temp_"+name+rand_sig+".fasta","w+")
	
	count_1,count_2,count_3=(0,0,0)
	f_count_2=0
	reads=dict()
	exact_threshold=num_of_seg
	over_threshold=num_of_seg+1
	under_threshold=num_of_seg-1
	for i in curr_result.iterkeys():
		if len(curr_result[i])>over_threshold:
			#print "ERROR_A",i, curr_result[i]
			f_count_2+=1
		elif len(curr_result[i])==over_threshold:
			#print "ERROR_C",i, curr_result[i]
			count_3+=1
		elif len(curr_result[i])==exact_threshold:
			reads[i]=curr_result[i][0][2]
		elif len(curr_result[i])==under_threshold:
			count_1+=1
		else:
			#print "ERROR_D",i, curr_result[i]
			count_2+=1

	print "ERROR_A:%d,ERROR_B:%d,ERROR_C:%d,ERROR_D:%d" % (f_count_2,count_1,count_3,count_2)
	del curr_result
	
	generate_fastas(source_file,reads,file_ptr)

	for i in file_ptr.iterkeys():
		file_ptr[i].close()

	BLAT=" "+options.local_align_path+" -out=pslx "+options.local_align_index
	for i in file_ptr.iterkeys():
		instrcution=BLAT+":"+i+" temp_"+i+rand_sig+".fasta  temp_"+i+rand_sig+".pslx "
		retcode = subprocess.call(instrcution, shell=True)
		if retcode < 0:
			print >>stderr, "Child was terminated by signal", -retcode
	
	print "Finish the running of BALT"

	tmp=argv[0].rfind("/")
	instrcution="cat  temp_*"+rand_sig+".pslx | tee "+output_file+".pslx |"+argv[0][:tmp+1]+"pslx2sam.pl  > "+output_file+".sam "
	retcode = subprocess.call(instrcution, shell=True)
	if retcode < 0:
		print >>stderr, "Child was terminated by signal", -retcode
	print "Finish coverting pslx format to SAM foramt."
	
	if debug==-1:
		instrcution="rm  temp_*"+rand_sig+"* "
		retcode = subprocess.call(instrcution, shell=True)
		if retcode < 0:
			print >>stderr, "Child was terminated by signal", -retcode
		print "Finsih cleaning temporary files."
	
	
