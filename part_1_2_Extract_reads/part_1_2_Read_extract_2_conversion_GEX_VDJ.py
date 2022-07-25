import os, sys, re
import glob
import multiprocessing
import pysam

# author: Juraj Michalik
# date: 22/07/15

# This demultiplexes files based on corrected barcodes
# uses pre-computation used by part 1, which is absolutely necessary for
# this to work properly

# class which stores info relevant to a given iteration:
# - run and lane files
# - iteration for which this is done
# - unique handles for input and output files that keep track of search
class Run_Lane:
	def __init__ (self, fastq_file_paths, fastq_target_file_paths, run, lane, iteration):
		
		# source files
		self.file_paths = fastq_file_paths
		self.files = []
		for path in self.file_paths: # load fastq files
			fastq_handle = open(path, 'r')
			self.files.append(fastq_handle)

		# target files
		self.target_file_paths = fastq_target_file_paths
		self.target_files = []
		for path in self.target_file_paths: # load fastq files
			target_fastq_handle = open(path, 'w')
			self.target_files.append(target_fastq_handle)
				
		self.run = run
		self.lane = lane
		self.iter = iteration

	# redefine print 
	def __str__(self):
		print_str = "Run_Lane_instance: \nFiles: " + ", ".join(self.file_paths) + ", ".join(self.target_file_paths)
		print_str += "\nRun: " + self.run + "\nLane: " + self.lane
		return(print_str)		

	# search for read in fastq (must be ordered before)
	# if both searched and searching ordered, we need to search only in one direction
	# done on all reads (R1, R2, possibly I1 and I2)
	def search_read_by_header(self, searched_header):
		while True: # ends when we hit end of file (should *not* happen or read is found)
			current_reads = []
			for self_file in self.files:
				actual_read = self_file.readline()
				if not actual_read:
					exit('Read header not found. This should not happen. Did you sort sam and/or fastqs?')
				else:
					for i in range(0,3): # three more lines
						actual_read += self_file.readline()
				
				current_reads.append(actual_read) 
			
			#print('@' + searched_header)
			#print(current_reads[0].split('\n')[0].split(' ')[0])
			#print('__')	
			
			# returns reads for all files for given run + lane combo, as list, each read is a string with whitespace conservation
			if ('@' + searched_header) == current_reads[0].split('\n')[0].split(' ')[0]:
				return(current_reads)

	# writes reads obtained by search function from Run_Lane_Source
	# must be in same order  !!!
	def write_reads(self, written_reads):
		for i, written_read in enumerate(written_reads):
			self.target_files[i].write(written_read)
		
			

# processes sam chunks and writes them to file by appartenance
def process_sam_chunks(sam_chunk, fastq_dirs_sorted, cell_barc_list, source_dict, result_dir, iter_n):
	
	# load sam
	samfile = pysam.AlignmentFile(sam_chunk, "r")
	
	# create target dir structure which will facilitate file writing
	# first create directory specific for this run (each parallel run will have its own directory for obvious reasons)
	iter_spec_result_dir = result_dir + os.sep + 'Fastq_iter_' + str(iter_n)
	
	# create directory where fastqs of this run only will be stored
	if not os.path.exists(iter_spec_result_dir):
		os.makedirs(iter_spec_result_dir)
		
	# now proceed with creation of target structure, using source structure
	target_dict = dict.fromkeys(list(source_dict.keys()))
	
	for run_key in target_dict:
		target_dict[run_key] = dict.fromkeys(list(source_dict[run_key].keys()))
		for lane_key in target_dict[run_key]:
			source_fastq_paths = source_dict[run_key][lane_key]
			
			# create run-specific sub-directories
			if not os.path.exists(iter_spec_result_dir + os.sep + run_key):
				os.makedirs(iter_spec_result_dir + os.sep + run_key)			
		
			target_fastq_paths = [iter_spec_result_dir + os.sep + run_key + os.sep + os.path.basename(x) for x in source_fastq_paths]
			target_dict[run_key][lane_key] = Run_Lane(source_fastq_paths, target_fastq_paths, run_key, lane_key, iter_n)
			#print(target_dict[run_key][lane_key])
	
	# now iterate over reads in samfile chunk - if they match desired gene (or are multi mapping) and barcode we keep them
	# if they are not mapping we keep them as well
	
	for sam_alignment in samfile.fetch():
		# search for tags
		# barcode first
		try: 
			feature_barcode_tag = sam_alignment.get_tag('fr')
			continue # means it's a feature barcoding read, so we skip
			
		except KeyError:
			pass
		
		
		try:
			cell_barcode = sam_alignment.get_tag('CB')

			if(cell_barcode in cell_barc_list):
									
				# we extract run ID and lane number from header and use them as keys for desired structures
				read_name_vals = sam_alignment.query_name.split(':')
				read_run = read_name_vals[2]
				read_lane = "L%03d"%int(read_name_vals[3])
					
				# now search for reads and write them to files
				matching_reads = target_dict[read_run][read_lane].search_read_by_header(sam_alignment.query_name)					
				target_dict[read_run][read_lane].write_reads(matching_reads)
												
		except KeyError:
			continue # bad barcodes are rejected (in the end, we demultiplex)
			
	samfile.close()

if __name__ == "__main__":
	sam_file_chunks_dir = sys.argv[1]			# sam chunks which will be analyzed
	fastq_dirs_sorted = sys.argv[2]				# sorted fastq dirs, can contian multiple split by comma
	cell_barc_list_path = sys.argv[3]			# list of cell barcodes
	result_dir = sys.argv[4]					# where results will be stored
	
	rep_path = os.path.dirname(os.path.normpath(sam_file_chunks_dir)) + os.sep + 'report_partII.txt'
	
	# writes parameter report to dir where bam is located
	target = open(rep_path, 'w')
	target.write('||USED PARAMETERS||\n')
	target.write('SAM CHUNKS: ' + sam_file_chunks_dir + '\n')
	target.write('SORTED FASTQS: ' + fastq_dirs_sorted + '\n')
	target.write('CELLULAR BARCODE LIST: ' + cell_barc_list_path + '\n')
	target.write('RESULT PATH: ' + result_dir + '\n')
	target.close()	
						
	# load barcode list as python object
	cell_barc_list_handle = open(cell_barc_list_path, 'r')
	cell_barc_list = cell_barc_list_handle.readlines()
	cell_barc_list_handle.close()
	
	# remove trailing whitespaces and quotes
	cell_barc_list = [x.strip().strip('"') + '-1' for x in cell_barc_list]
	
	# load source fastq file in dictionary - this will be handy later for each 
	# prep and creation of keys of each level of dir
	fastq_dirs_sorted = fastq_dirs_sorted.split(',')
	fastq_dir_keys = [os.path.basename(os.path.normpath(x)).split('_')[0] for x in fastq_dirs_sorted]
	
	source_dict = dict.fromkeys(fastq_dir_keys)
	
	# prepare dictionary itself with Run_Lane_Source objects at the end
	for i, fastq_dir in enumerate(fastq_dirs_sorted):
		fastq_files = sorted(glob.glob(fastq_dir + os.sep + '*'))
		
		lane_keys = []
		
		# detect how many lanes are there
		for fastq_file in fastq_files:
			lane_key = re.search('_L(\d)+_', fastq_file).group(0).strip('_')
			if lane_key not in lane_keys:
				lane_keys.append(lane_key)
				
		source_dict[fastq_dir_keys[i]] = dict.fromkeys(lane_keys)
		
		# now create Run_Lane_Source object for each lane
		for lane_key in source_dict[fastq_dir_keys[i]]:
			fastq_files_by_lane = [x for x in fastq_files if lane_key in x]
			source_dict[fastq_dir_keys[i]][lane_key] = fastq_files_by_lane	
	
	# create directory where finalized fastqs will be stored
	if not os.path.exists(result_dir):
		os.makedirs(result_dir)
	
	# run core of this program in parallel manner
	sam_file_chunks = sorted(glob.glob(sam_file_chunks_dir + os.sep + '*'))
	
	proc_list = []
	
	for iter_n,sam_chunk in enumerate(sam_file_chunks):
		print('Started iteration ' + str(iter_n) + '.')
		proc_run = multiprocessing.Process(target = process_sam_chunks, args = (sam_chunk, fastq_dirs_sorted, cell_barc_list, source_dict, result_dir, iter_n,))
		
		proc_list.append(proc_run)
		
	for proc in proc_list:
		proc.start()
	
	for proc in proc_list: # wait until everything is completed
		proc.join()	
		
	# merge resulting files
	# create directory with merged results
	merged_fastq_path = result_dir + os.sep + 'Merged'
	if not os.path.exists(merged_fastq_path):
		os.makedirs(merged_fastq_path)
		
	for run_key in source_dict:
		read_types = ['I1', 'I2', 'R1', 'R2']
		for lane_key in source_dict[run_key]:
			for read_type in read_types:
				glob_path = result_dir + os.sep + 'Fastq_iter*' + os.sep + run_key + os.sep
				iter_res_files = sorted(glob.glob(glob_path + '*' + lane_key + '_' + read_type + '*'))
				
				merge_cmd_line = ' '
				if len(iter_res_files) > 0:
					for iter_res_file in iter_res_files:
						
						merge_cmd_line += iter_res_file
						merge_cmd_line += ' '
		
					merged_file_basename = os.path.basename(iter_res_files[0])	
					merge_cmd_line = 'cat' + merge_cmd_line + ' > ' + result_dir + os.sep + 'Merged' 
					merge_cmd_line += os.sep + run_key + '_' + merged_file_basename
			
					print(merge_cmd_line)
					os.system(merge_cmd_line)
		

	
