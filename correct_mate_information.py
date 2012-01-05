import argparse



def bundle_process(bundle,output_file_handler,validator):
	ID=bundle[0]
	if len(bundle[1])==1:
		bundle[1][0][1]=str(int(bundle[1][0][1]) | 9) # 1+8
		output_file_handler.write("%s\n" % '\t'.join(bundle[1][0]))
	elif len(bundle[1])==2:
		rule1= (bundle[1][0][2]==bundle[1][1][2]) # Chr. == Chr.
		rule2= abs(int(bundle[1][0][3])-int(bundle[1][1][3]))<=50000 # default value of Tophat
		rule3= validator(int(bundle[1][0][1]),int(bundle[1][1][1]))  # augmented validator for seqnecing protocol
		
		# If read_id contains mate index information, i.e.,
		#	AAAA_1,AAAA_2, or
		#	XXXX/1,XXXX/2
		# Then we can identified the order of mates.		
		if bundle[1][0][0][-2:] < bundle[1][1][0][-2:]:
			first=bundle[1][0]
			second=bundle[1][1]
		else:
			first=bundle[1][1]
			second=bundle[1][0]
		
		
		if True==(rule1 and rule2 and rule3):	
			#The alignment of this paired-end read is mated
			if bundle[1][0][0][-2:] ==  bundle[1][1][0][-2:]:
				#the index of the segment in the template is unknown. unset both 0x40 and 0x80
				second[1] = str(int(second[1]) | 3) # 2+1
				first[1] = str(int(first[1]) | 3)  # 2+1
			else:
				second[1] = (int(second[1]) | 131) # 128+2+1
				if second[1] & 16 == 16:
					first[1] = str(int(first[1]) | 99 ) # 64+32+2+1
				else:
					first[1] = str(int(first[1]) | 67)  # 64+2+1
				second[1] = str(second[1])
				
				first[0]=first[0][:-2]
				second[0]=second[0][:-2]
		else:
			if bundle[1][0][0][-2:] ==  bundle[1][1][0][-2:]:
				#the index of the segment in the template is unknown. unset both 0x40 and 0x80
				second[1] = str(int(second[1]) | 1) # 1
				first[1] = str(int(first[1]) | 1)  # 1
			else:
				second[1] = (int(second[1]) | 129) # 128+1
				if second[1] & 16 == 16:
					first[1] = str(int(first[1]) | 97 ) # 64+32+1
				else:
					first[1] = str(int(first[1]) | 65)  # 64+1
				second[1] = str(second[1])
		
		output_file_handler.write("%s\n" % '\t'.join(first))
		output_file_handler.write("%s\n" % '\t'.join(second))
				
	else:
		return False	
	return True
	

def correct_mate_information(input_file,output_file,mate_information):
	fr=open(input_file,'r')
	fw=open(output_file,'w+')
	if mate_information=="--fr":
		is_mated=lambda x,y: ((x ^ y) & 16) == 16
	elif mate_information=="--ff":
		is_mated=lambda x,y: ((x ^ y) & 16) == 0
	
	bundle=('',[])
	for line in fr:
		line=line.strip("\r\n").split("\t")
		if (int(line[1])& 4)==4:
			continue
		ID=line[0][:-2]
		if len(bundle[1])==0:
			bundle=(ID,[])
		elif ID!=bundle[0]:
			bundle_process(bundle,fw,is_mated)
			bundle=(ID,[])
		bundle[1].append(line)
		
	
	bundle_process(bundle,fw,is_mated)
	fw.close()
	fr.close()
	return True


if __name__ == '__main__':
	description='''correct the mate information of SAM format'''
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', dest='input_file', metavar='PATH', type=str, help='',default='')
	parser.add_argument('-o', dest='output_file', metavar='PATH', type=str, help='',default='')
	parser.add_argument('-m', dest='mate_info', metavar='STR', type=str, nargs='*',help='',default='--fr')
	
	args = parser.parse_args()
	
	correct_mate_information(args.input_file,args.output_file,args.mate_info)
	
	
