#!/usr/bin/env python

'''
RNASEQRp is a version using pipe to connect different counterparts of RNASEQR. 

'''

import argparse
import subprocess
from time import clock, gmtime, strftime,time
from os import devnull
import sys
from selected_the_multiply_mapped_reads_of_the_same_genomic_position import *
import bundle_paired_end as bundle_pe
import correct_mate_information as cmi
#import SPEEDUP

def bowtie_unique(input_file):
	tmp=input_file.split(".")
	tmp.append("unique."+tmp.pop())
	output_file_u=".".join(tmp)
	fw_u=open(output_file_u,'w+')
	
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
	f.close()
	return True	


def load_parameters_from_config(arguments):
	'''
	This function load the parameters into the arguments object.
	Since the RNASEQR would give commandline parameter higher priority than configuration file, the function would modify the arguments object if the the value of the the argument is identical to default value.
	'''
	f=open(arguments.configuration,'r')
	for line in f:
		line=line.strip('\t\r\n ')
		if line=='' or line[0]=='#':
			continue
		
		line=line.split('=')
		if line[0] in arguments:
			tmp=getattr(arguments,line[0])
			if tmp=='' or tmp == -1:
				setattr(arguments,line[0],line[1])
	f.close()
	return True


if __name__ == '__main__':
	description='''RNASEQR is a nucleotide sequence mapper/aligner and is designed specifically for RNA-seq data analysis. It takes advantage of the annotated transcripts and genomic reference sequences to obtain high quality mapping/alignment results. For any inquiry, please contact Leslie Chen PhD (lchen AT systemsbiology DOT org)'''
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('configuration', metavar='setting.config', type=str, nargs='?',help='configuration file. Please note that the command-line arguments will override the parameters specified in this file.',default='')
	parser.add_argument('-a', dest='aligner_path', metavar='PATH', type=str, help='aligner_path',default='')
	parser.add_argument('-m', dest='align_mismatch', metavar='Integer', type=int, help='align_mismatch',default=-1)
	parser.add_argument('-b', dest='seg_align_mismatch', metavar='Integer', type=int, help='segment align_mismatch',default=-1)
	parser.add_argument('-p', dest='number_of_thread', metavar='Integer', type=int, help='number_of_thread',default=-1)
	parser.add_argument('-s', dest='phred_type', metavar='TYPE', type=str, help='phred_type: "Sanger","Solexa","Illumina1.3","Illumina1.5+"',choices=["Sanger","Solexa","Illumina1.3","Illumina1.5+"],default='Sanger')
	parser.add_argument('-l', dest='local_aligner_path', metavar='PATH', type=str, help='local_aligner_path',default='')
	parser.add_argument('-n', dest='genome_index_prefix_for_local_aligner', metavar='PATH', type=str, help='genome_index_prefix_for_local_aligner',default='')
	parser.add_argument('-t', dest='annotation_type', metavar='TYPE', type=str, help='annotation_type',choices=["BED","GTF"],default='')
	parser.add_argument('-i', dest='annotation_index_prefix', metavar='PATH', type=str, help='annotation_index_prefix',default='')
	parser.add_argument('-q', dest='input_fastq', metavar='PATH', type=str, help='Input FASTQ file. If the input is paired-end data, use comma to separate them in order.',default='')
	parser.add_argument('-f', dest='annotation_file', metavar='PATH', type=str, help='annotation_file',default='')
	parser.add_argument('-g', dest='genome_index_prefix', metavar='PATH', type=str, help='genome_index_prefix',default='')
	parser.add_argument('-r', dest='short_read_length', metavar='Integer', type=int, help='The length of short read',default=-1)
	parser.add_argument('-e', dest='num_of_hit_anchor_to_be_valid', metavar='Integer', type=int, help='The number of hit anchors to be valid',default=-1)
	parser.add_argument('-x', dest='anchor_length', metavar='Integer', type=int, help='The length of anchors',default=-1)
	parser.add_argument('-C', dest='colorspace', metavar='Integer', type=int, help='If use colorspace, set the value of this option as 1.',default=-1)
	parser.add_argument('-G', dest='paired_end_gap', metavar='String', type=str, help='[Paired-end Only] The insertion length of paired-end reads. The format is INT-INT, e.g., 200-400.',default="")
	parser.add_argument('-o', dest='max_num_overlap_exon', metavar='Integer', type=str, help='Maximum number of overlaping exons in annotation',default=40)
	parser.add_argument('--version', action='version', version='RNASEQR 1.0, 2011-12-01.')
	parser.add_argument('--debug', dest='debug', metavar='Integer', type=int, help='Set value 1 to open the debug mode',default=-1)
	args = parser.parse_args()
	
	#initialize parameters
	RNASEQR_path,aligner_name,aligner_path,phred_type,annotation_type,annotation_file='','','','','',''
	local_aligner_path,genome_index_prefix_for_local_aligner='',''
	annotation_index_prefix,input_fastq,genome_index_prefix='','',''
	short_read_length,align_mismatch,seg_align_mismatch,number_of_thread = -1,-1,-1,-1
	colorspace,paired_end_gap,debug,num_of_hit_anchor_to_be_valid = -1,-1,-1,-1
	anchor_length,max_num_overlap_exon=-1,-1
	paired_end_input=False
	
	if args.configuration!='':
		#load configuration
		if True != load_parameters_from_config(args):
			print >>sys.stderr, 'Fail to load the confiuration file.'
	
	#load augments
	if RNASEQR_path =='':
		tmp = sys.argv[0].rfind("/")
		if tmp != -1:
			RNASEQR_path = sys.argv[0][:tmp]
		else:
			RNASEQR_path = sys.argv[0]
	
	if aligner_name =='':
		aligner_name="Bowtie"
	
	aligner_path=args.aligner_path
	local_aligner_path=args.local_aligner_path
	align_mismatch=int(args.align_mismatch)
	seg_align_mismatch=int(args.seg_align_mismatch)
	number_of_thread=int(args.number_of_thread)
	short_read_length=int(args.short_read_length)	
	annotation_type=args.annotation_type
	annotation_file=args.annotation_file
	annotation_index_prefix=args.annotation_index_prefix
	input_fastq=args.input_fastq
	genome_index_prefix=args.genome_index_prefix
	genome_index_prefix_for_local_aligner=args.genome_index_prefix_for_local_aligner
	colorspace=int(args.colorspace)
	paired_end_gap=args.paired_end_gap
	debug=int(args.debug)
	num_of_hit_anchor_to_be_valid=int(args.num_of_hit_anchor_to_be_valid)
	anchor_length=int(args.anchor_length)
	max_num_overlap_exon=int(args.max_num_overlap_exon)
	
	if annotation_type=='':
		annotation_type='GTF'
	
	# Splicing short read to generate anchors.
	if anchor_length==-1:
		# By default, it will be adjusted by RNASEQR
		if short_read_length >= 75:
			anchor_length=25
		elif short_read_length <75:
			anchor_length=16
	if num_of_hit_anchor_to_be_valid==-1:
		num_of_hit_anchor_to_be_valid=(short_read_length/anchor_length)-1
	
	debug_logger_backup=sys.stderr
	if debug==1:
		sys.stderr=devnull
	elif debug==2:
		sys.stderr=open('RNASEQR.log','w+')
	
	
	phred_type = args.phred_type
	bowtie_phred={"Sanger":"","Solexa":"--solexa-quals","Illumina1.3":"--solexa1.3-quals","Illumina1.5+":"--solexa1.3-quals"}
	if phred_type == '':
		phred_type="Sanger"
	Phred = phred_type
	phred_type = bowtie_phred[phred_type]
	
	
	#show_the_parameters
	print "\n=== Sepicified Parameters ==="
	print "RNASEQR path: %s" % RNASEQR_path
	print "Aligner: %s [%s]" % (aligner_name,aligner_path)
	print "\tMismatch: %d\n\tSegment mismatch:%d\n\tThread number:%d" % (align_mismatch,seg_align_mismatch,number_of_thread)
	print "\tPhred: %s\n\tAnnotation: %s [%s,%s]\n\tGenome: [%s]" % (Phred,annotation_type,annotation_file,annotation_index_prefix,genome_index_prefix)
	print "\tMaximum number of overlapping exons in annotation: [%s]" % (max_num_overlap_exon)
	
	print "Local aligner: %s " % local_aligner_path
	print "Local aligner index: %s " % genome_index_prefix_for_local_aligner
	print "Input: %s" % input_fastq
	print "\tRead Length: %d" % short_read_length
	
	if num_of_hit_anchor_to_be_valid<=0:
		# No reasonable num_of_hit_anchor_to_be_valid
		print "Anchor-and-Map: Disabled"
	else:
		print "Anchor-and-Map: Enabled"
	print "\tAnchor Length: %d" % anchor_length
	print "\tnum_of_hit_anchor_to_be_valid: %d" % num_of_hit_anchor_to_be_valid
	
	
	if colorspace==-1:
		print "Colorspace: False"
	else:
		print "Colorspace: True"
	print "==============================\n\n"
	
	
	# ====================================================================
	#
	# Main procedure
	#
	# ====================================================================
	if input_fastq.find(",")!=-1:
		paired_end_input=True
		input_fastq=input_fastq.split(",")
		paired_end_gap=map(lambda x:int(x),paired_end_gap.split("-"))
	
	if colorspace==-1:
		colorspace=""
	else:
		colorspace=" -C --col-keepends "
	
	print "Starting Preporcessing: ",strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
	ENST_structures=extract_gene_structure(annotation_file,annotation_type)
	print "Finish Preporcessing: ",strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
	
	
	# ================================
	# Stage 1: Algin to Transcriptome
	# ================================
	print "Start Aligning: ",strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
	
	if paired_end_input==True:
		mate_orientation="--fr"
		if colorspace!="":
			mate_orientation="--ff"
		
		instrcution='''%s  -v  %d  -p  %d %s -k %d -m %d --best --strata -S --sam-nohead -I %d -X %d --un phase1.pe.unmapped.fastq --max /dev/null %s %s -S --sam-nohead  %s  -1 %s -2 %s ''' % (aligner_path,align_mismatch,number_of_thread,mate_orientation,max_num_overlap_exon,max_num_overlap_exon,paired_end_gap[0],paired_end_gap[1],phred_type,colorspace,annotation_index_prefix,input_fastq[0],input_fastq[1])
		print instrcution
		prev_time=time()
		bundle_pe.verify_for_paired_end_reads(instrcution,"phase1.unique.sam",ENST_structures)
	else:
		instrcution='''%s  -v  %d  -p  %d  -k %d -m %d -q --un phase1.unmapped.fastq --max /dev/null --best --strata  %s %s -S --sam-nohead  %s  %s  ''' % (aligner_path,align_mismatch,number_of_thread,max_num_overlap_exon,max_num_overlap_exon,phred_type,colorspace,annotation_index_prefix,input_fastq)
		prev_time=time()
		#SPEEDUP.verification_for_bowtie(instrcution,"phase1.unique.sam",ENST_structures)
		verification_for_bowtie_PIPE(instrcution,"phase1.unique.sam",ENST_structures)
	
	print "Phase 1 alignment time cost: ", (time() - prev_time), "sec."
	
	
	# ================================
	# Stage 2: Algin to genome
	# ===============================

	if paired_end_input==True:
		#Bowtie output all unmappable paired-end reads into one file.
		
		instrcution="cat phase1.pe.unmapped* > phase1.unmapped.fastq "
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "WRONG: ",instrcution
		
		instrcution='''%s  -v  %d  -p  %d  -k 1 -m 1 -q --un phase2.unmapped.fastq --max /dev/null --best --strata  %s %s -S --sam-nohead  %s  %s  > phase2.sam ''' % (aligner_path,align_mismatch,number_of_thread,phred_type,colorspace,genome_index_prefix,"phase1.unmapped.fastq")
		prev_time=time()
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "phase 2 error"
			print instrcution
		print "Phase 2 alignment time cost: ", (time() - prev_time), "sec."
	else:
		instrcution='''%s  -v  %d  -p  %d  -k 1 -m 1 -q --un phase2.unmapped.fastq --max /dev/null --best --strata  %s %s -S --sam-nohead  %s  %s  > phase2.sam ''' % (aligner_path,align_mismatch,number_of_thread,phred_type,colorspace,genome_index_prefix,"phase1.unmapped.fastq")
		print instrcution
		prev_time=time()
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "phase 2 error"
			print instrcution
		print "Phase 2 alignment time cost: ", (time() - prev_time), "sec."
	
	
	# ================================
	# Stage 3: local alignment
	# ================================
	print "num_of_hit_anchor_to_be_valid",num_of_hit_anchor_to_be_valid
	if num_of_hit_anchor_to_be_valid >= 0:
		print "num_of_hit_anchor_to_be_valid",num_of_hit_anchor_to_be_valid
		'''
		RNASEQR can not find reasonable num_of_hit_anchor_to_be_valid. 
		Thus we just skip the Stage 3.
		'''
		instrcution='''python '''+RNASEQR_path+'/generate_anchors.py -a '+str(anchor_length)+' -r '+str(short_read_length)+' phase2.sam phase2.anchors.fastq'
		print instrcution
		prev_time=time()
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "generate stage 2 input reads"
			print instrcution
		print "generate stage 2 input reads time cost: ", (time() - prev_time), "sec."
		
		instrcution='''%s  -v  %d  -p  %d  -k %d -m %d -q --best --strata  %s  -S --sam-nohead  %s ''' % (aligner_path,seg_align_mismatch,number_of_thread,max_num_overlap_exon,max_num_overlap_exon,annotation_index_prefix,"phase2.anchors.fastq")
		prev_time=time()
		verification_for_bowtie_PIPE(instrcution,"stage2.T.unique.sam",ENST_structures)
		
		instrcution='''%s  -v  %d  -p  %d  -k 1 -m 1 -q --best --strata  %s  -S --sam-nohead  %s > stage2.G.sam ''' % (aligner_path,seg_align_mismatch,number_of_thread,genome_index_prefix,"phase2.anchors.fastq")
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "stage 2 alignemnt (G)"
			print instrcution
		print "Stage 2 alignment time cost: ", (time() - prev_time), "sec."

		print "Finish Aligning: ",strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		print "Start Processing Alignment Result",strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

		prev_time=time()
		bowtie_unique("phase2.sam")
		print "Extract Uniqueness (1) Time cost: ", (time() - prev_time), "sec."
		
		prev_time=time()
		bowtie_unique("stage2.G.sam")
		print "Extract Uniqueness (2) Time cost: ", (time() - prev_time), "sec."
		
		print "Start Alginement Refinment",strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		
		instrcution=' python '+RNASEQR_path+'/extract_sequence_for_stage2.py phase2.unmapped.fastq stage2.T.unique.sam stage2.G.unique.sam stage2.report -d '+str(debug)+' -n '+str(num_of_hit_anchor_to_be_valid)+' -i '+genome_index_prefix_for_local_aligner+' -p '+local_aligner_path
		print instrcution
		prev_time=time()
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "Anchored and Align"
			print instrcution
		print "Anchored and Align Time cost: ", (time() - prev_time), "sec."

		#compute the allowed mismatch for spliced alignment
		allowed_trim_size_of_spliced_alignment=str((short_read_length-2)-1)
		instrcution=' python '+RNASEQR_path+'/analyze_gaps_of_alignments.py -f -g 0 -m '+allowed_trim_size_of_spliced_alignment+' stage2.report.sam >  stage2.report.m'+allowed_trim_size_of_spliced_alignment+'g0.sam'
		prev_time=time()
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "Extract Uniqueness (1)"
			print instrcution
		print "Extract Uniqueness (1) Time cost: ", (time() - prev_time), "sec."
		
		prev_time=time()
		bowtie_unique("stage2.report.m"+allowed_trim_size_of_spliced_alignment+"g0.sam")
		print "Extract Uniqueness (2) Time cost: ", (time() - prev_time), "sec."

		prev_time=time()
		instrcution=' python '+RNASEQR_path+'/convert_blat_resulted_sam_to_standard_sam.py stage2.report.m'+allowed_trim_size_of_spliced_alignment+'g0.unique.sam phase2.unmapped.fastq > final.sam '
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "Extract Uniqueness (3)"
			print instrcution
		print "Convert BLAT Time cost: ", (time() - prev_time), "sec."
	
	
	
	# ================================
	# Stage 4: Preparing output and clean up temperory files
	# ================================
	if paired_end_input==True:
		instrcution=' cat final.sam phase2.unique.sam | sort -k 1,1 > stage4.sam '
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "Extract Uniqueness (3)"
		# correcting the mate information
		cmi.correct_mate_information('stage4.sam','final.sam',mate_orientation)
		instrcution=' cat phase1.unique.sam >> final.sam '
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "Extract Uniqueness (3)"
		print instrcution
	else:
		instrcution=' cat phase1.unique.sam phase2.unique.sam >> final.sam '
		prev_time=time()
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "Extract Uniqueness (3)"
		print instrcution
		print "Catencation Time cost: ", (time() - prev_time), "sec."

	if debug==1:
		sys.stderr=debug_logger_backup
	elif debug==2:
		sys.stderr.close()
		sys.stderr=debug_logger_backup
	else:
		instrcution=' rm -f stage* stage* phase* '
		prev_time=time()
		retcode = subprocess.call(instrcution,shell=True)
		if retcode < 0:
			print >>sys.stderr, "remove temp. files: stage*  stage* phase* "
			print instrcution
	

