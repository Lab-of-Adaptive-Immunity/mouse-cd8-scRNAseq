import sys, os, glob
import pysam
# @author: Juraj Michalik
# @date: 02/01/2021
# performs VDJ search on GEX data and returns found VDJ table
#
# @VDJ CONFIGURATION:
# - you must have mixcr software and precise its path (see main of this script)
# (https://github.com/milaboratory/mixcr/releases - we used version 3.0.12)	
# - you have to dowload bamtofastq program from 10X site
# (https://cf.10xgenomics.com/misc/bamtofastq-1.2.0)
# 
# these files should:
# - be filtered to genes of interest only (otherwise this will take forever and a day)
# - must be indexed

# reads and splits bam file to 
def split_bam(bam_path, resdir_path, minread = 25):
	sam = pysam.AlignmentFile(bam_path, 'rb')
	sam_iter = sam.fetch()
	
	# step one: split in lists
	barcode_dict = {}
	splitted_bam_list = []
	barcode_counter = 0
	for alignment in sam_iter:
		barcode = alignment.get_tag('CR')
		if barcode in barcode_dict:
			splitted_bam_list[barcode_dict[barcode]].append(alignment)
		
		else:
			barcode_dict[barcode] = barcode_counter
			splitted_bam_list.append([alignment])
			barcode_counter += 1
			
	# step two: write all in files, one file opened at a time
	# ignore all cells that have less than minread reads per cell: it will be impossible to get any info from low-read cells
	for cell in splitted_bam_list:
		if len(cell) > minread:
			writeh = pysam.AlignmentFile(resdir_path + "/bam/bam_" + cell[0].get_tag('CR') + ".bam", "wb", template=sam)
			for alignment in cell:
				writeh.write(alignment)
			
			writeh.close()

# transforms all files 		
def bamtofastq(resdir_path, bamtofastq10x_path):
	print('Now converting to fastq')
	filelist = glob.glob(resdir_path+'/bam/*')
	for bamfile in filelist:
		root = bamfile.split('/')[-1].split('.')[0].split('_')[-1]
		dest_path = resdir_path+'/fastq/'+root+'/'
		os.system(bamtofastq10x_path + ' ' + bamfile + ' '+dest_path)
		
		# fusion all R1 and R2 files resulting from bamtofastq
		# R1:
		R1s = sorted(glob.glob(dest_path + '/*/*R1*'))
		if(len(R1s) > 0):
			os.system('cat '+' '.join(R1s)+' > '+dest_path+'/fastq_R1.fastq.gz')
		
		# R2:
		R2s = sorted(glob.glob(dest_path + '/*/*R2*'))
		if(len(R2s) > 0):
			os.system('cat '+' '.join(R2s)+' > '+dest_path+'/fastq_R2.fastq.gz')	
		print(R1s)
		print(R2s)

# runs mixcr on all split bam files
def run_mixcr(resdir_path, mixcr_path, an_species):
	fastqs = glob.glob(resdir_path+'/fastq/*')	
	for fastq_dir in fastqs:
		reads_1 = fastq_dir + '/' + 'fastq_R1.fastq.gz'
		reads_2 = fastq_dir + '/' + 'fastq_R2.fastq.gz'
		root = fastq_dir.split('/')[-1].split('.')[0].split('_')[-1]
		target = resdir_path + '/mixcr/' + root
		if not os.path.exists(target):
			os.makedirs(target)
		os.system(mixcr_path + ' analyze shotgun -s ' + an_species + ' --starting-material rna --only-productive --receptor-type tcr '
		 + reads_1 + ' ' + reads_2 + ' ' + target + '/' + root)  
		
				
# combine mixcr results of all cells in one global output
def generate_final_output(resdir_path):
	mixcr_results = glob.glob(resdir_path+'/mixcr/*')
	# read header from first file, but add 
	headerfile = open(mixcr_results[0]+'/'+ mixcr_results[0].split('/')[-1] + '.clonotypes.TRA.txt', 'r')
	header = headerfile.readline()
	headerfile.close()
	header = 'barcode\t'+header
	
	trac = open(resdir_path+'/GEX_TRA.csv', 'w')
	trbc = open(resdir_path+'/GEX_TRB.csv', 'w')
	trac.write(header)
	trbc.write(header)
	
	for result_dir in mixcr_results:
		barcode = result_dir.split('/')[-1]
		# read TCRA and TCRB results
		
		if os.path.exists(result_dir+'/'+barcode+'.clonotypes.TRA.txt'):
			TRAh = open(result_dir+'/'+barcode+'.clonotypes.TRA.txt', 'r')
			TRA_l = TRAh.readlines()
			TRAh.close()
			
			for line in TRA_l[1:]:
				trac.write(barcode+'\t'+line)
		
		if os.path.exists(result_dir+'/'+barcode+'.clonotypes.TRB.txt'):
			TRBh = open(result_dir+'/'+barcode+'.clonotypes.TRB.txt', 'r')
			TRB_l = TRBh.readlines()
			TRBh.close()
			
			for line in TRB_l[1:]:
				trbc.write(barcode+'\t'+line)
		
	trac.close()
	trbc.close()


if __name__ == "__main__":	
	# config
	bamtofastq10x_path = 'bamtofastq-1.2.0'
	mixcr_path = 'mixcr-3.0.12/mixcr'
	
	# params
	bam_path = sys.argv[1]
	resdir_path = sys.argv[2]
	an_species = 'mmu'
	if(len(sys.argv) == 4):
		an_species = sys.argv[3]
	
	if not os.path.exists(resdir_path + '/bam'):
		os.makedirs(resdir_path + '/bam')
		
	if not os.path.exists(resdir_path + '/fastq'):
		os.makedirs(resdir_path + '/fastq')

	if not os.path.exists(resdir_path + '/mixcr'):
		os.makedirs(resdir_path + '/mixcr')
	
	split_bam(bam_path, resdir_path)
	bamtofastq(resdir_path, bamtofastq10x_path)
	run_mixcr(resdir_path, mixcr_path, an_species)
	generate_final_output(resdir_path)
