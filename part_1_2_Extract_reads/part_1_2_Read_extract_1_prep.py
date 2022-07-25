import os, sys
import glob

# author: Juraj Michalik
# date: 22/07/15

# This prepares bam and corresponding fastqs for demultiplexing by barcodes
# in following way:
# - sorts fastq	and bam files by header name (to speed up the process)
# - splits bam to multiple sub-bams to allow parallel processing

# requires imported python packages and samtools

if __name__ == "__main__":
	bam_file = sys.argv[1]						# where bam is stored
	bam_file_chunks_dir = sys.argv[2]			# dir where chunks of sorted bam will be stored
	fastq_dirs = sys.argv[3]					# dir with fastqs, multiple dirs can be entered (separated by comma)
	fastq_dirs_sorted = sys.argv[4]				# dir where sorted fastqs will be put, multiple dirs can be entered separated by comma and same len. as above (will create if does not exist)
	chunk_len = sys.argv[5]						# length of single chunk to which bam file will be split
	ncores = sys.argv[6]						# number of cores to use for samtools view
	
	rep_path = os.path.dirname(bam_file) + os.sep + 'report.txt'
	
	# writes parameter report to dir where bam is located
	target = open(rep_path, 'w')
	target.write('||USED PARAMETERS||\n')
	target.write('BAM: ' + bam_file + '\n')
	target.write('BAM CHUNKS: ' + bam_file_chunks_dir + '\n')
	target.write('FASTQS SOURCE: ' + fastq_dirs + '\n')
	target.write('FASTQS TARGET: ' + fastq_dirs_sorted + '\n')
	target.write('BAM CHUNKS LENGTH: ' + chunk_len + '\n')
	target.write('SAMTOOLS NCORES: ' + ncores + '\n')
	target.close()
	
	
	# sort bam by read names
	bam_name_sort = os.path.splitext(bam_file)[0] + '_sort.bam'
	os.system('samtools sort -@ ' + ncores + ' -n ' + bam_file + ' > '+ bam_name_sort)
	
	fastq_dirs = fastq_dirs.split(',')
	fastq_dirs_sorted = fastq_dirs_sorted.split(',')
	
	if(len(fastq_dirs) != len(fastq_dirs_sorted)):
		sys.exit("Source and result fastq directory lists are of inequal length!")
	
	for i in range(len(fastq_dirs)):
	
		# sort fastqs by name
		fastq_paths = glob.glob(fastq_dirs[i] + os.sep + '*')
		if not os.path.exists(fastq_dirs_sorted[i]):
			os.makedirs(fastq_dirs_sorted[i])
	
		for fastq_path in fastq_paths:
			fastq_path_sorted =  fastq_dirs_sorted[i] + os.sep + os.path.splitext(os.path.basename(fastq_path))[0]
			print(fastq_path)
			print(fastq_path_sorted)
			os.system('zcat ' + fastq_path + ' | paste - - - - | sort -k1,1 -V --parallel ' + ncores + ' -t " " | tr "\t" "\n" > ' + fastq_path_sorted)

	# split bam to multiple sub-files by chunk size
	# uses samtools view for conversion to sam then split

	if not os.path.exists(bam_file_chunks_dir):
		os.makedirs(bam_file_chunks_dir)

	sam_name = os.path.splitext(bam_name_sort)[0] + '.sam'
	split_prefix = bam_file_chunks_dir + os.sep + os.path.splitext(os.path.basename(bam_name_sort))[0]
	print(sam_name)
	print(split_prefix)
	os.system('samtools view -@ ' + ncores + ' ' + bam_name_sort + ' > ' + sam_name) # convert to sam
	os.system('split -d -l ' + chunk_len + ' ' + sam_name + ' ' + split_prefix) # split to chunks
	
	# re-add header
	# get header for sorted bam
	os.system('samtools view -H ' + bam_name_sort + ' > bam_header')
	sam_chunks = glob.glob(bam_file_chunks_dir + os.sep + '*')
	
	# add header
	for chunk in sam_chunks:
		os.system('cat bam_header ' + chunk + ' > ' + chunk + 'tmp')
		os.system('mv ' + chunk + 'tmp ' + chunk) # move at the place of original file

	os.system('rm bam_header')
