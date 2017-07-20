#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import glob
import urllib2
import re
import gzip
import os
from subprocess import call
from StringIO import StringIO


NCBI_GENBANK_ASSEMBLY_URL = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
OUTPUT_PATH = '/s/angus/index/databases/kraken_databases/ncbi_genbank_assembly_summary.txt'
OUTPUT_FASTA_FOLDER = '18July2017_kraken_db'
FTP_REGEXP = re.compile(r'(reference genome|representative genome)')
PLASMID_REGEXP = re.compile(r'plasmid')


class KrakenParser(object):
	def __init__(self):
		self.ncbi_data = {}
		self.ftp_addresses = {}
		self.temp_folder = '/'.join(OUTPUT_PATH.split('/')[:-1]) + '/temporary_file_dir'
		self.outfolder = '/'.join(OUTPUT_PATH.split('/')[:-1]) + '/' + OUTPUT_FASTA_FOLDER
		self.total_files = 0

	def parse_ncbi_data(self, data_path):
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
			#print (extension)
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
	f = urllib2.urlopen(url)
	data = f.read()
	with open(output_path, 'wb') as out:
		out.write(data)


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
