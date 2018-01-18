#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import glob
import urllib2
import re
import gzip
import os
from subprocess import call
from StringIO import StringIO

# The following variables in uppercase must be edited for your server and to specify which genomes you want to include.
OUTPUT_PATH = '/s/angus/index/databases/kraken_databases/combined_assembly_summary.txt'
OUTPUT_FASTA_FOLDER = 'meglab-customkraken-Jan2017'
GROUP_list = ["bacteria","archaea","viral"]
#archaea|bacteria|fungi|invertebrate|metagenomes|other|plant|protozoa|vertebrate_mammalian|vertebrate_other|viral 	# These are all the options for genome categories. 

# These variables can be edited as desired
NCBI_GENBANK_ASSEMBLY_URL = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/'
FTP_REGEXP = re.compile(r'(reference genome|representative genome)')
PLASMID_REGEXP = re.compile(r'plasmid')

class KrakenParser(object):
	def __init__(self):
		self.ncbi_data = {}
		self.ftp_addresses = {}
		self.temp_folder = '/'.join(OUTPUT_PATH.split('/')[:-1]) + '/temporary_file_dir'
		self.outfolder = '/'.join(OUTPUT_PATH.split('/')[:-1]) + '/' + OUTPUT_FASTA_FOLDER
		self.total_files = 0

	def parse_ncbi_data(self, data_path): # This takes in the assembly_summary files
		with open(data_path, 'rb') as f:
			data = f.read().split('\n')
			for entry in data:
				if not entry:
					continue
				if entry[0] == '#':
					continue
				line = entry.split('\t')
				self.ncbi_data.setdefault(line[0], line[1:])

	def find_reliable_ftp_addresses(self):
		for acc, vals in self.ncbi_data.items():
			if FTP_REGEXP.match(vals[3]):
				self.ftp_addresses.setdefault(acc, (vals[18], vals[4]))
	
	def download_reliable_genomes(self):
		if not os.path.exists(self.temp_folder):
			os.makedirs(self.temp_folder)
		total_files = len(self.ftp_addresses.keys())
		print_progress_bar(0, total_files, prefix='Progress: ', suffix='Complete', length=60)
		counter = 0
		plasmids = 0
		for acc, vals in self.ftp_addresses.items():
			extension = vals[0].split('/')[-1]
			counter += 1
			request = urllib2.Request(vals[0] + '/' + extension + '_genomic.fna.gz')
			request.add_header('Accept-encoding', 'gzip')
			response = urllib2.urlopen(request)
			buf = StringIO(response.read())
			f = gzip.GzipFile(fileobj=buf)
			data = {k: v for k, v in fasta_parse(f)}
			with open(self.temp_folder + '/' + acc + '.fasta', 'wb') as out:
				for header, seq in data.items():
					if PLASMID_REGEXP.match(header):
						plasmids += 1
						out.write('>kraken:taxid|32630|{}|{}\n{}\n'.format(acc, header, seq))
					else:
						out.write('>kraken:taxid|{}|{}|{}\n{}\n'.format(vals[1], acc, header, seq))
			print_progress_bar(counter, total_files, prefix='Progress: ', suffix='Complete', length=60)
		print('\n')
		print('{} plasmid sequences downloaded out of {} total genomes'.format(plasmids, counter))
		self.total_files = counter
	
	def run_dust(self):
		print_progress_bar(0, self.total_files, prefix='Progress: ', suffix='Complete', length=60)
		if not os.path.exists(self.outfolder):
			os.makedirs(self.outfolder)
		counter = 0
		for f in glob.glob(self.temp_folder + '/*.fasta'):
			counter += 1
			filename = f.split('/')[-1]
			call(['dustmasker', '-in', f, '-out', self.outfolder + '/' + filename, '-infmt', 'fasta', '-outfmt', 'fasta'])
			print_progress_bar(counter, self.total_files, prefix='Progress: ', suffix='Complete', length=60)
			print('\n')

	def mask_low_complexity_regions(self):
		print_progress_bar(0, self.total_files, prefix='Progress: ', suffix='Complete', length=60)
		counter = 0
		for f in glob.glob(self.outfolder + '/*.fasta'):
			counter += 1
			with open(f, 'r+') as fastafile:
				data = {k: v for k, v in fasta_parse(fastafile)}
				fastafile.seek(0)
				for header, seq in data.items():
					newstring = ''
					for nuc in seq:
						if nuc.islower():
							newstring += 'N'
						else:
							newstring += nuc
					fastafile.write('>{}\n{}\n'.format(header, newstring))
				fastafile.truncate()
			print_progress_bar(counter, self.total_files, prefix='Progress: ', suffix='Complete', length=60)
			print('\n')
			

def download_genbank_file(url, output_path):
	group_data = []
	for i in GROUP_list:
		url_i = (url + i + "/assembly_summary.txt")
		f = urllib2.urlopen(url_i)
		print(url_i)
		data = f.read()
		group_data.append(data)
	with open(output_path, 'wb') as out:
		for i in group_data:
			out.write(i)

def fasta_parse(fasta_file):
	# Skip whitespace
	while True:
		line = fasta_file.readline()
		if line is "":
			return  # Empty file or premature end of file?
		if line[0] is ">":
			break
	while True:
		if line[0] is not ">":
			raise ValueError("Records in FASTA should begin with '>'")
		header = line[1:].rstrip()
		all_lines = []
		line = fasta_file.readline()
		while True:
			if not line:
				break
			if line[0] is ">":
				break
			all_lines.append(line.rstrip())
			line = fasta_file.readline()
		yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
		if not line:
			return  # Stop Iteration
	assert False, "Should not reach this line"


def print_progress_bar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100):
	"""
	Call in a loop to create terminal progress bar
	@params:
		iteration   - Required  : current iteration (Int)
 		total       - Required  : total iterations (Int)
		prefix      - Optional  : prefix string (Str)
		suffix      - Optional  : suffix string (Str)
		decimals    - Optional  : positive number of decimals in percent complete (Int)
		length      - Optional  : character length of bar (Int)
		fill        - Optional  : bar fill character (Str)
	"""
	percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
	filledLength = int(length * iteration // total)
	bar = '█' * filledLength + '-' * (length - filledLength)
	sys.stdout.write('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix))
	# Print New Line on Complete
	if iteration == total: 
		print()
		
		
if __name__ == '__main__':
	if len(sys.argv) > 1:
		if sys.argv[1] == 'new':
			download_genbank_file(NCBI_GENBANK_ASSEMBLY_URL, OUTPUT_PATH)
	print('Parsing NCBI data file...\n')
	kp = KrakenParser()
	kp.parse_ncbi_data(OUTPUT_PATH)
	kp.find_reliable_ftp_addresses()

	print('Downloading genomes...\n')
	kp.download_reliable_genomes()

	print('Running DUST...\n')
	kp.run_dust()

	print('Masking low complexity regions...\n')
	kp.mask_low_complexity_regions()

	print('Done.\n')

## Work in progress section, these tasks are currently done manually and should be implemented into the code
	
## # download UniVec database to the library, and modify the headers to be picked up by Centrifuge and Kraken
# Thanks to email collaboration with Steven Salzberg and Florian Breitwieser
#wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec | sed 's/>/>kraken:taxid|32630| /' > UniVec.fa
# download EmVec database
#wget ftp://ftp.ebi.ac.uk/pub/databases/emvec/emvec.dat.gz
#gunzip -c emvec.dat | sed -n '/^DE/,/^\/\//p' | grep -E '^DE|^    ' | sed '/^  /s/  *[0-9]*$//g;;/^  /s/ //g;s/^DE/>kraken:taxid|32630| /;/^[^>]/s/.*/\U&/' > emvec.fa
#Furthermore (and quite importantly), get the Enterobacteria phage PhiX-174 genome from http://www.ncbi.nlm.nih.gov/nuccore/CP004084.1 and also modify the header to start with '>kraken:taxid|32630| '. The PhiX-174 in the vector database is a bit different from that one, for some reason.
#cat UniVec.fa emvec.fa phiX174 > vector_contaminants.fasta
#move that vector_contaminants.fasta into the OUTPUT_FASTA_FOLDER prior to adding the  genomes to the database, and before "building" the database

## Use code to add all genomes to new database location, then download taxonomy and finally build the database

#for file in temp_folder/*.fasta
#do
    #kraken-build --add-to-library $file --db $DBNAME
#done
#kraken-build --download-taxonomy --db meglab_kraken_db_Jan2018
#kraken-build --build --db meglab_kraken_db_Jan2018

## Need to erase the temporary folder and run the clean task for kraken
kraken-build --db $DBNAME --clean

## Do we need to " re-run a kraken script after the database building, to set the lowest common ancestor (LCA) for all k-mers, that appear in contaminant sequences, to the artificial vector taxid (32630). For this I had to patch some files in kraken, and I can share this with you. Normally, those k-mers that are seen in both contaminants and genomes get a taxid 1 (i.e. root)  - which makes it more difficult to filter artificial reads.?"

