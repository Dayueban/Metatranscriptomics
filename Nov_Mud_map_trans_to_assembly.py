#!/usr/bin/python

#PYTHON SCRIPT 
#written by: Richard Wolfe
#
#to run type: python single_copy_genes.py -i <inputfile> -o <outputfile> -m <min sequence length>
#         or: ./single_copy_genes.py attributes
#
#   if error: /usr/bin/python^M: bad interpreter: No such file or directory
#      -there is a windows endl after shebang
#      -open in vi 
#         once in vi type:
#           :set ff=unix<return>
#           :x<return>
#
#
#   reads all the sam file in the directory and make a table
# 
#


import sys      #for exit command and maxint
import argparse #to get command line args 
                #needed to install argparse module because using python 2.6
                #and argparse comes with python 2.7
                #  sudo easy_install argparse
import os       #to run system commands
#import datetime #to make timestamp
import glob     #for * in file name
import toolbox


def make_raw_reads_table():

	#Note: multiple alignments are counted as a read that maps

	dir_files = []
	dir_files.append('Plant_11_14_A')
	dir_files.append('Plant_11_14_B')
	dir_files.append('Plant_11_14_C')
	dir_files.append('Mud_11_14_A')
	dir_files.append('Mud_11_14_B')
	dir_files.append('Mud_11_14_C')
	dir_files.append('Plant_9_15_A')
	dir_files.append('Plant_9_15_B')
	dir_files.append('Plant_9_15_C')
	dir_files.append('Mud_9_15_A')
	dir_files.append('Mud_9_15_B')
	dir_files.append('Mud_9_15_C')

	reads = [] #list of lists
	ids = []

	for item in dir_files:
		file_list = glob.glob('mismatches_*_' + item + '*.sam')   #put all  files in the list

		if len(file_list) == 0 or len(file_list) > 1:
			print "Error ... Not only 1 file"
			sys.exit(1)

		print file_list[0]

		header_lines = 0
		data_lines = 0
		mapped_reads = 0

		temp = []
		for i in ids:
			temp.append(0)

		f = open(file_list[0], "rU")
		line = f.readline()
		while line:
			if line.startswith( '@' ):
       				header_lines += 1
			else: #this is a data line
				data_lines += 1
		
				cols = line.split()
				flags = int(cols[1])  #get flags

				if not flags & 4:  #if aligns
					mapped_reads += 1
				
					scaff = cols[2]
					
					#print line
					#print scaff
					#sys.exit(0)

					if scaff in ids:
						index = ids.index(scaff)
						temp[index] += 1
					else:
						ids.append(scaff)  #add this id
		
						#add 0 fpkm to all other samples
						for j in reads:
							j.append(0)

						#add 1 to this list
						temp.append(1)

					


			line = f.readline()

		f.close()

		reads.append(temp)

		print "header lines = ", header_lines
		print "data lines = ", data_lines
		print "mapped reads = ",mapped_reads
		print ""

		#break #temp so only 1 sample



	#write to table
	args.make_raw_reads_table.write("gene_id")
	for item in dir_files:
		args.make_raw_reads_table.write("\t" + item)
	 
	args.make_raw_reads_table.write("\n")

	i = 0
	while i < len(ids):
		args.make_raw_reads_table.write(ids[i])
		for item in reads:
			args.make_raw_reads_table.write("\t" + str(item[i]))
		args.make_raw_reads_table.write("\n")

		i += 1


	print ""
	print "Script finished"
	sys.exit(0)
	

def make_table():

	#need the fasta file that was used to run bowtie
	#if not args.make_table_fasta:
	#	print "Error ... need --make_table_fasta fasta file used to make bowtie index"
	#	sys.exit(1)

	if not args.folder:
		print "Error ... need --folder extension after sample name"
		sys.exit(1)

	#get the name of the output file
	out_file_name = ""

	if args.make_table:
		args.make_table.close()
		out_file_name = args.make_table.name
	elif args.make_rounded_table:
		args.make_rounded_table.close()
		out_file_name = args.make_rounded_table.name
	else:
		#should not get here
		print "Error ... no output file for make_table"
		sys.exit()

	out_file = open(out_file_name, "w")
	

	#make a list of all the ids in the fasta file
	ids = []
	samples = []
	fpkm_list = []  #this is a list of lists

	#line = args.make_table_fasta.readline()
	#while line:
	#	if line.startswith(">"):
	#		line = line.rstrip()
	#		line = line[1:] #remove >
	#		#print line
	#		ids.append(line.split()[0]) #get 1st space
	#	line = args.make_table_fasta.readline()

	#args.make_table_fasta.close()

	#print "Number of ids in fasta file = ",len(ids)
			

	#dir_files = glob.glob('/ORG-Data/Wetlands/Metatranscripts/*')   #put all  files in the list
	dir_files = []
	dir_files.append('Plant_11_14_A')
	dir_files.append('Plant_11_14_B')
	dir_files.append('Plant_11_14_C')
	dir_files.append('Mud_11_14_A')
	dir_files.append('Mud_11_14_B')
	dir_files.append('Mud_11_14_C')
	dir_files.append('Plant_9_15_A')
	dir_files.append('Plant_9_15_B')
	dir_files.append('Plant_9_15_C')
	dir_files.append('Mud_9_15_A')
	dir_files.append('Mud_9_15_B')
	dir_files.append('Mud_9_15_C')
	
	for item in dir_files:
        	#there are 3 Plant 11 samples, 3 plant 9, 3Mud 11, 3 Mud ...
		#print item   #Ex: /ORG-Data/Wetlands/Metatranscripts/Plant_11_14_A

		#extract the sample name
		sample_name = item.split("/")[-1]
		#print sample_name
		samples.append(sample_name)

		folder_name = sample_name + args.folder
		print folder_name

		if not os.path.isdir(folder_name):
			print "Error ... Path is not a directory"
			sys.exit(1)

		#make a temp list of fpkm values
		temp = []
		for i in ids:
			temp.append("0")

		f = open(folder_name + "/genes.fpkm_tracking", "rU")

		line = f.readline() #first line is header
		line = f.readline()

		if not line: #file is empty
			print "File is empty:  " + folder_name

		while line:
			line = line.rstrip()
			cols = line.split("\t")

			gene_id = cols[6]
			#there is a :1-86.. on end we need to remove
			gene_id = gene_id.split(":")[0]
			
			fpkm = cols[9]

			if args.make_rounded_table:  #if we are rounding the FPKM value
				i = int(round(float(fpkm)))  #makes a float from the string, rounds it and then make int
				fpkm = str(i)  #make it a string

			if gene_id in ids:
				index = ids.index(gene_id)
				temp[index] = fpkm
			else:
				ids.append(gene_id)  #add this id
		
				#add 0 fpkm to all other samples
				for j in fpkm_list:
					j.append("0")

				#add 0 to this list
				temp.append(fpkm)
			

			#print gene_id
			#print fpkm
			#sys.exit(0)

			line = f.readline()

		f.close()

		fpkm_list.append(temp)

	#write to table
	out_file.write("gene_id")
	for item in samples:
		out_file.write("\t" + item)
	 
	out_file.write("\n")

	i = 0
	while i < len(ids):
		out_file.write(ids[i])
		for item in fpkm_list:
			out_file.write("\t" + item[i])
		out_file.write("\n")

		i += 1
	

	print ""
	print "Number of samples = ", len(samples)
	print ""
	print "Script finished"
	sys.exit(0)



def run_database():
	#close the fasta file we only want the name of the file
	args.fasta.close()

	file_name = args.fasta.name.split("/")[-1]  #if path get filename

	print "file_name = ",file_name

	#make a bowtie database with the fasta file
	if args.skip_bowtie == "F":
		print "Making bowtie index"
		cmd = "bowtie2-build " + args.fasta.name + " " + file_name
		toolbox.run_system(cmd)


	#dir_files = glob.glob('/ORG-Data/Wetlands/Metatranscripts/*')   #put all  sam files in the list
	dir_files = []
	dir_files.append('Plant_11_14_A')
	dir_files.append('Plant_11_14_B')
	dir_files.append('Plant_11_14_C')
	dir_files.append('Mud_11_14_A')
	dir_files.append('Mud_11_14_B')
	dir_files.append('Mud_11_14_C')
	dir_files.append('Plant_9_15_A')
	dir_files.append('Plant_9_15_B')
	dir_files.append('Plant_9_15_C')
	dir_files.append('Mud_9_15_A')
	dir_files.append('Mud_9_15_B')
	dir_files.append('Mud_9_15_C')
	
	reads = []  #number of reads in sam file
	mapped_reads = [] #number of mapped reads in sam file


	#print the  files
	for sample_name in dir_files:
        	#there are 3 Plant 11 samples, 3 plant 9, 3Mud 11, 3 Mud ...
		print sample_name   #Ex: /ORG-Data/Wetlands/Metatranscripts/Plant_11_14_A

		#extract the sample name
		#sample_name = item.split("/")[-1]
		#sample_name = item
		#print sample_name

		#if not sample_name.startswith("Mud_11"): #only want to process Nov Mud
		#	continue
	
		#if sample_name.startswith("Mud_11"): #want to process all exceptNov Mud
		#	continue
		
		#see if already processed
		#dont_do_list = ["Mud_11_14_C","Mud_9_15_A","Mud_9_15_C","Plant_11_14_A","Plant_11_14_B","Plant_9_15_A","Plant_9_15_B"]

		#if sample_name in dont_do_list:
      		#	continue

		bt_db = file_name 

		sam_file = sample_name + "_mappedto_" + bt_db + ".sam"
		#r1 = "../" + sample_name + "/R1_no_Ribo_trimmed.fastq"
		#r2 = "../" + sample_name + "/R2_no_Ribo_trimmed.fastq"
		r1 = "/home2/projects/Wetlands/Metatranscripts/" + sample_name + "/R1_no_Ribo_trimmed.fastq"
		r2 = "/home2/projects/Wetlands/Metatranscripts/" + sample_name + "/R2_no_Ribo_trimmed.fastq"
		mis_match_file = "mismatches_" + str(args.mismatches) + "_" + sam_file
		bam_file = mis_match_file + ".bam"
		sorted_bam_file = "SORTED_" + bam_file 
		cufflinks_dir = sample_name + "_cufflinks_" + bt_db + "_mis_" + str(args.mismatches)
		corrected_cufflinks_dir =  sample_name + "_cufflinks_corrected_" + bt_db + "_mis_" + str(args.mismatches)

		#bowtie to assembly with multiple align -a option
		if args.skip_bowtie == "F":
			cmd = "bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -a -p 20 -x " + bt_db + " -S " + sam_file + " -1 " + r1 + " -2 " + r2
			toolbox.run_system(cmd)

		#change reads with mismatches <= 2
		cmd = "python /ORG-Data/scripts/sam_file.py -i " + sam_file + " -v " + str(args.mismatches) + " -o " + mis_match_file
		toolbox.run_system(cmd)

		#convert to bam
		cmd = "samtools view -@ 20 -bS " + mis_match_file + " > " + bam_file
		toolbox.run_system(cmd)

		#sort bam
		cmd = "samtools sort -@ 20 " + bam_file + " " + sorted_bam_file
		toolbox.run_system(cmd)
	
		#NOTE extra .bam???
		#runn cufflinks
		sorted_bam_file = sorted_bam_file + ".bam"
		cmd = "/home2/opt/Cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks -o " + cufflinks_dir + " " + sorted_bam_file
		toolbox.run_system(cmd)

		#runn cufflinks with corrected for multialign
		cmd = "/home2/opt/Cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks -u -o " + corrected_cufflinks_dir + " " + sorted_bam_file
		toolbox.run_system(cmd)


	
	print ""
	print "Script finished"
	sys.exit(0)

def run_qc():

	dir_files = glob.glob('/ORG-Data/Wetlands/Metatranscripts/*')   #put all  sam files in the list
	reads = []  #number of reads in sam file
	mapped_reads = [] #number of mapped reads in sam file


	#print the  files
	for item in dir_files:
        	#there are 3 Plant 11 samples, 3 plant 9, 3Mud 11, 3 Mud ...
		print item   #Ex: /ORG-Data/Wetlands/Metatranscripts/Plant_11_14_A

		#extract the sample name
		sample_name = item.split("/")[-1]
		print sample_name

		#if sample_name.startswith("Mud_11"): #only want to process if not Nov Mud
		#	continue

			
		#make a dir with sample name
		#if directory already made then it will go to next entry
		try:
			os.mkdir(sample_name)
		except:
			continue


		#unzip reads with ribo removed
		#gunzip -c /ORG-Data/Wetlands/Metatranscripts/10376.5.159127.TGACCA.anqrpht.fastq.gz > 10376.5.159127.TGACCA.anqrpht.fastq
		unzipped_file = sample_name + "/" + sample_name + ".anqrpht.fastq"	
		cmd = "gunzip -c " + item + "/*.anqrpht.fastq.gz > " + unzipped_file
		toolbox.run_system(cmd)
		#unzipped file is Plant_11_14_A/Plant_11_14_A.anqrpht.fastq
	

		#seperate reads into R1_no_Ribo.fastq and R2_no_Ribo.fastq
		#-these are interleaved forward and reverse reads
		#paste - - - - - - - - < Plant_11_14_A.anqrpht.fastq  | tee >(cut -f 1-4 | tr "\t" "\n" > R1_no_Ribo.fastq)  | cut -f 5-8 | tr "\t" "\n" > R2_no_Ribo.fastq
		#cmd = "paste - - - - - - - - < " + unzipped_file + ' | tee >(cut -f 1-4 | tr "\\t" "\\n" > ' + sample_name + '/R1_no_Ribo.fastq) | cut -f 5-8 | tr "\\t" "\\n" > ' + sample_name + "/R2_no_Ribo.fastq"
		#toolbox.run_system(cmd)
		toolbox.deinterleave_fastq_reads(unzipped_file, sample_name + "/R1_no_Ribo.fastq", sample_name + "/R2_no_Ribo.fastq")
	
		#trim the reads
		cmd = "sickle pe -f " + sample_name + "/R1_no_Ribo.fastq -r " + sample_name + "/R2_no_Ribo.fastq -t sanger -o " + sample_name + "/R1_no_Ribo_trimmed.fastq -p " + sample_name + "/R2_no_Ribo_trimmed.fastq -s " + sample_name + "/R1R2_no_Ribo_trimmed.fastq"
		toolbox.run_system(cmd)
		
	
	print ""
	print "Script finished"
	sys.exit(0)



def run_jordan_mud():

	dir_files = glob.glob('/ORG-Data/Wetlands/Metatranscripts/*')   #put all  sam files in the list
	reads = []  #number of reads in sam file
	mapped_reads = [] #number of mapped reads in sam file


	#print the  files
	for item in dir_files:
        	#there are 3 Plant 11 samples, 3 plant 9, 3Mud 11, 3 Mud ...
		print item   #Ex: /ORG-Data/Wetlands/Metatranscripts/Plant_11_14_A

		#extract the sample name
		sample_name = item.split("/")[-1]
		print sample_name

		#if not sample_name.startswith("Mud_11"): #only want to process Nov Mud
		#	continue
	
		if sample_name.startswith("Mud_11"): #want to process all exceptNov Mud
			continue

		bt_db = "Transcript_Mapping_DatabaseInputTest_Jordan"

	

		#bowtie to assembly with multiple align -a option
		cmd = "bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -a -p 40 -x " + bt_db + " -S " + sample_name + "_mappedto_DatabaseInputTest.sam -1 ../" + sample_name + "/R1_no_Ribo_trimmed.fastq -2 ../" + sample_name + "/R2_no_Ribo_trimmed.fastq"
		toolbox.run_system(cmd)

		#change reads with mismatches <= 2
		cmd = "python /ORG-Data/scripts/sam_file.py -i " + sample_name + "_mappedto_DatabaseInputTest.sam -v 2 -o " + sample_name + "_mismatches_2_mappedto_DatabaseInputTest.sam"
		toolbox.run_system(cmd)

		#convert to bam
		cmd = "samtools view -@ 60 -bS " + sample_name + "_mismatches_2_mappedto_DatabaseInputTest.sam > " + sample_name + "_mismatches_2_mappedto_DatabaseInputTest.bam"
		toolbox.run_system(cmd)

		#sort bam
		cmd = "samtools sort -@ 60 " + sample_name + "_mismatches_2_mappedto_DatabaseInputTest.bam " + sample_name + "_mismatches_2_mappedto_DatabaseInputTest_SORTED.bam"
		toolbox.run_system(cmd)
	
		#NOTE extra .bam???
		#runn cufflinks
		cmd = "/home2/opt/Cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks -o " + sample_name + "cufflinks_DatabaseInputTest " + sample_name + "_mismatches_2_mappedto_DatabaseInputTest_SORTED.bam.bam"
		toolbox.run_system(cmd)

		#runn cufflinks with corrected for multialign
		cmd = "/home2/opt/Cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks -u -o " + sample_name + "cufflinks_corrected_DatabaseInputTest " + sample_name + "_mismatches_2_mappedto_DatabaseInputTest_SORTED.bam.bam"
		toolbox.run_system(cmd)

	print ""
	print "Script finished"
	sys.exit(0)

def run_garrett_mud():

	dir_files = glob.glob('/ORG-Data/Wetlands/Metatranscripts/*')   #put all  sam files in the list
	reads = []  #number of reads in sam file
	mapped_reads = [] #number of mapped reads in sam file


	#print the  files
	for item in dir_files:
        	#there are 3 Plant 11 samples, 3 plant 9, 3Mud 11, 3 Mud ...
		print item   #Ex: /ORG-Data/Wetlands/Metatranscripts/Plant_11_14_A

		#extract the sample name
		sample_name = item.split("/")[-1]
		print sample_name

		#if not sample_name.startswith("Mud_11"): #only want to process Nov Mud
		#	continue

		if sample_name.startswith("Mud_11"): # want to process all except Nov Mud
			continue
	
		bt_db = "NovMethanotrophBin"

	

		#bowtie to assembly with multiple align -a option
		cmd = "bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -a -p 40 -x " + bt_db + " -S " + sample_name + "_mappedto_NovMethanotrophBin.sam -1 ../" + sample_name + "/R1_no_Ribo_trimmed.fastq -2 ../" + sample_name + "/R2_no_Ribo_trimmed.fastq"
		toolbox.run_system(cmd)

		#change reads with mismatches <= 2
		cmd = "python /ORG-Data/scripts/sam_file.py -i " + sample_name + "_mappedto_NovMethanotrophBin.sam -v 2 -o " + sample_name + "_mismatches_2_mappedto_NovMethanotrophBin.sam"
		toolbox.run_system(cmd)

		#convert to bam
		cmd = "samtools view -@ 60 -bS " + sample_name + "_mismatches_2_mappedto_NovMethanotrophBin.sam > " + sample_name + "_mismatches_2_mappedto_NovMethanotrophBin.bam"
		toolbox.run_system(cmd)

		#sort bam
		cmd = "samtools sort -@ 60 " + sample_name + "_mismatches_2_mappedto_NovMethanotrophBin.bam " + sample_name + "_mismatches_2_mappedto_NovMethanotrophBin_SORTED.bam"
		toolbox.run_system(cmd)
	
		#NOTE extra .bam???
		#runn cufflinks
		cmd = "/home2/opt/Cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks -o " + sample_name + "cufflinks_NovMethanotrophBin " + sample_name + "_mismatches_2_mappedto_NovMethanotrophBin_SORTED.bam.bam"
		toolbox.run_system(cmd)

		#runn cufflinks with corrected for multialign
		cmd = "/home2/opt/Cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks -u -o " + sample_name + "cufflinks_corrected_NovMethanotrophBin " + sample_name + "_mismatches_2_mappedto_NovMethanotrophBin_SORTED.bam.bam"
		toolbox.run_system(cmd)

	print ""
	print "Script finished"
	sys.exit(0)


print "Script started ..."


#create an argument parser object
#description will be printed when help is used
parser = argparse.ArgumentParser(description='A script to run AMPHORA2')

#add the available arguments -h and --help are added by default
#if the input file does not exist then program will exit
#if output file does not exit it will be created
# args.input is the input file Note: cant write to this file because read only
# args.output is the output file
# args.m is the minimum seq length
#parser.add_argument('-o', '--output_file', type=argparse.FileType('w'), help='output file', required=True)
#parser.add_argument('-i', '--id', help='id', required=True)
parser.add_argument('-q', '--qc', help='run qc on all files except Mud_11_14  T or F', default="F")
parser.add_argument('--jm', help='run jordans db with mud_11  T or F', default="F")
parser.add_argument('--gm', help='run garretts db with mud_11  T or F', default="F")
parser.add_argument('-f','--fasta',type=argparse.FileType('rU'), help='make a bowtie index from this file and run all reads through cufflinks')
parser.add_argument('-t','--make_table',type=argparse.FileType('w'),help='make a table from all the cufflinks folders in this directory')
parser.add_argument('--make_rounded_table',type=argparse.FileType('w'),help='make a table with rounded FPKM values from all the cufflinks folders in this directory')
parser.add_argument('--make_raw_reads_table',type=argparse.FileType('w'),help='make a table with number of raw reads from all the mismatches .sam files in this directory')

#parser.add_argument('--make_table_fasta',type=argparse.FileType('rU'), help='fasta file that was used to make the bowtie index')
parser.add_argument('--folder',help='folder extension after the sample name used to make a table')
parser.add_argument('-v','--mismatches',type=int,help='number of mismatches') #required=True

parser.add_argument('--skip_bowtie',help='skip bowtie step  T or F',default="F") #required=True

#get the args
args = parser.parse_args()

#additional argument tests

#print "length of args = ",len(sys.argv)
#print sys.argv
#sys.exit(0)

#check if arguments
if len(sys.argv) == 1: #the first argument is the script name
	print "Error .. No arguments were supplied"
	sys.exit(1)

if args.make_table and args.make_rounded_table:
	print "Error ... can't have both args --make_table and --make_rounded_table"
	sys.exit(1)

if args.make_table == None and args.make_rounded_table == None and args.make_raw_reads_table == None:  #if not making table then mismatches is reqd and must be > 0
	if args.mismatches == None:
		print "Error ... -v mismatches is required"
		sys.exit(1)
	elif args.mismatches < 0:
		print "Error .. -v is less than 0"
		sys.exit(1)

if args.make_raw_reads_table:
	make_raw_reads_table()

if args.fasta: #if the fasta argument is supplied
	if args.skip_bowtie != "T":
		if args.skip_bowtie != "F":
			print "--skip_bowtie must be T or F"
			sys.exit(1)
	run_database()

if args.qc != "T":
	if args.qc != "F":
		print "ERROR ... -q must be T or F"
		sys.exit(1)
if args.qc == "T":
	run_qc()

if args.jm != "T":
	if args.jm != "F":
		print "ERROR ... --jm must be T or F"
		sys.exit(1)
if args.jm == "T":
	run_jordan_mud()

if args.gm != "T":
	if args.gm != "F":
		print "ERROR ... --jm must be T or F"
		sys.exit(1)
if args.gm == "T":
	run_garrett_mud()


if args.make_table or args.make_rounded_table:
	make_table()

#Test print the args
#print args



#check if files exist from a previous time
#if os.path.isfile('*.aln') or os.path.isfile('*.mask') or os.path.isfile('*.pep') or os.path.isfile('*.phylotype'):
#	print("Files exist from previous run. Please delete or move to a different folder")
#        sys.exit(0)

dir_files = glob.glob('/ORG-Data/Wetlands/Metatranscripts/*')   #put all  sam files in the list
reads = []  #number of reads in sam file
mapped_reads = [] #number of mapped reads in sam file


#if len(sam_files) == 0:
#	print("Error .. There are no sam files in this dir")
#	sys.exit(0)


#print the  files
for item in dir_files:
        #there are 3 Plant 11 samples, 3 plant 9, 3Mud 11, 3 Mud ...
	print item   #Ex: /ORG-Data/Wetlands/Metatranscripts/Plant_11_14_A

	#extract the sample name
	sample_name = item.split("/")[-1]
	print sample_name
	
	if not sample_name.startswith("Mud_11"): #only want to process Nov Mud
		continue

	
	#make a dir with sample name
	os.mkdir(sample_name)

	
	#unzip reads with ribo removed
	#gunzip -c /ORG-Data/Wetlands/Metatranscripts/10376.5.159127.TGACCA.anqrpht.fastq.gz > 10376.5.159127.TGACCA.anqrpht.fastq
	unzipped_file = sample_name + "/" + sample_name + ".anqrpht.fastq"	
	cmd = "gunzip -c " + item + "/*.anqrpht.fastq.gz > " + unzipped_file
	toolbox.run_system(cmd)
	#unzipped file is Plant_11_14_A/Plant_11_14_A.anqrpht.fastq
	

	#seperate reads into R1_no_Ribo.fastq and R2_no_Ribo.fastq
	#-these are interleaved forward and reverse reads
	#paste - - - - - - - - < Plant_11_14_A.anqrpht.fastq  | tee >(cut -f 1-4 | tr "\t" "\n" > R1_no_Ribo.fastq)  | cut -f 5-8 | tr "\t" "\n" > R2_no_Ribo.fastq
	#cmd = "paste - - - - - - - - < " + unzipped_file + ' | tee >(cut -f 1-4 | tr "\\t" "\\n" > ' + sample_name + '/R1_no_Ribo.fastq) | cut -f 5-8 | tr "\\t" "\\n" > ' + sample_name + "/R2_no_Ribo.fastq"
	#toolbox.run_system(cmd)
	toolbox.deinterleave_fastq_reads(unzipped_file, sample_name + "/R1_no_Ribo.fastq", sample_name + "/R2_no_Ribo.fastq")
	
	#trim the reads
	cmd = "sickle pe -f " + sample_name + "/R1_no_Ribo.fastq -r " + sample_name + "/R2_no_Ribo.fastq -t sanger -o " + sample_name + "/R1_no_Ribo_trimmed.fastq -p " + sample_name + "/R2_no_Ribo_trimmed.fastq -s " + sample_name + "/R1R2_no_Ribo_trimmed.fastq"
	toolbox.run_system(cmd)


	#map the trimmed reads to the assembly db with multiple mappings
	bt_db = ""

	if sample_name.startswith("Plant_11"):
		bt_db = "bowtie_databases/Nov_plant_scaffold.fa"
	elif sample_name.startswith("Plant_9"):
		bt_db = "bowtie_databases/Aug_plant_scaffold.fa"
	elif sample_name.startswith("Mud_11"):
		bt_db = "bowtie_databases/Nov_mud_scaffold.fa"
	elif sample_name.startswith("Mud_9"):
		bt_db = "bowtie_databases/Aug_mud_scaffold.fa"
	else:
		print "Bowtie index error"
		sys.exit(1)

	#bowtie to assembly with multiple align -a option
	cmd = "bowtie2 --fast -a -p 40 -x " + bt_db + " -S " + sample_name + "/mappedto_assembly.sam -1 " + sample_name + "/R1_no_Ribo_trimmed.fastq -2 " + sample_name + "/R2_no_Ribo_trimmed.fastq"
	toolbox.run_system(cmd)

	#run to first mcrA_database_ALL.txt db
	bt_db = "bowtie_databases/mcrA_database_ALL.txt"
	cmd = "bowtie2 --fast -a -p 40 -x " + bt_db + " -S " + sample_name + "/mcrA_database_ALL.sam -1 " + sample_name + "/R1_no_Ribo_trimmed.fastq -2 " + sample_name + "/R2_no_Ribo_trimmed.fastq"
	toolbox.run_system(cmd) 

	#run to first methanogen scaffolds  db
	bt_db = "bowtie_databases/Nov_Mud_CombinedScaffolds_forDB.txt"
	cmd = "bowtie2 --fast -a -p 40 -x " + bt_db + " -S " + sample_name + "/Nov_Mud_CombinedScaffolds_forDB.sam -1 " + sample_name + "/R1_no_Ribo_trimmed.fastq -2 " + sample_name + "/R2_no_Ribo_trimmed.fastq"
	toolbox.run_system(cmd) 
	
		
	#sys.exit(0)  #stop after 1
	

#combine all the reads and map to the 3 bowtie databases
sample_name = "combined_Mud_11"
os.mkdir(sample_name)

#cat the trimmed reads
cmd = "cat Mud_11_14_A/R1_no_Ribo_trimmed.fastq Mud_11_14_B/R1_no_Ribo_trimmed.fastq Mud_11_14_C/R1_no_Ribo_trimmed.fastq > combined_Mud_11/R1_no_Ribo_trimmed.fastq"
toolbox.run_system(cmd)

cmd = "cat Mud_11_14_A/R2_no_Ribo_trimmed.fastq Mud_11_14_B/R2_no_Ribo_trimmed.fastq Mud_11_14_C/R2_no_Ribo_trimmed.fastq > combined_Mud_11/R2_no_Ribo_trimmed.fastq"
toolbox.run_system(cmd)

#bowtie to assembly with multiple align -a option
bt_db = "bowtie_databases/Nov_mud_scaffold.fa"
cmd = "bowtie2 --fast -a -p 40 -x " + bt_db + " -S " + sample_name + "/mappedto_assembly.sam -1 " + sample_name + "/R1_no_Ribo_trimmed.fastq -2 " + sample_name + "/R2_no_Ribo_trimmed.fastq"
toolbox.run_system(cmd)

#run to first mcrA_database_ALL.txt db
bt_db = "bowtie_databases/mcrA_database_ALL.txt"
cmd = "bowtie2 --fast -a -p 40 -x " + bt_db + " -S " + sample_name + "/mcrA_database_ALL.sam -1 " + sample_name + "/R1_no_Ribo_trimmed.fastq -2 " + sample_name + "/R2_no_Ribo_trimmed.fastq"
toolbox.run_system(cmd) 

#run to first methanogen scaffolds  db
bt_db = "bowtie_databases/Nov_Mud_CombinedScaffolds_forDB.txt"
cmd = "bowtie2 --fast -a -p 40 -x " + bt_db + " -S " + sample_name + "/Nov_Mud_CombinedScaffolds_forDB.sam -1 " + sample_name + "/R1_no_Ribo_trimmed.fastq -2 " + sample_name + "/R2_no_Ribo_trimmed.fastq"
toolbox.run_system(cmd) 

print ""
print "Number of  files in directory = ", len(dir_files)

print ""
print "Script finished..."
