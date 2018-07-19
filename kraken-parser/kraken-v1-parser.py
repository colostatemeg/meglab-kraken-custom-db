## Kraken version 1 kraken parser to remove plasmids
import sys
import re

PLASMID_REGEXP = re.compile(r'plasmid')

def fasta_parse(infile):
    with open(infile, 'r') as fastaFile:
        # Skip whitespace
        while True:
            line = fastaFile.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            allLines = []
            line = fastaFile.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                allLines.append(line.rstrip())
                line = fastaFile.readline()
            yield header, "".join(allLines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration
        assert False, "Should not reach this line"
		
def rename_plasmid_seqs(inseqs, outseqs):
		data = {k: v for k, v in fasta_parse(inseqs)}
		with open(outseqs , 'w') as out:
			plasmids = 0
			non_plasmids = 0
			for header, seq in data.items():
				print('This is working!')
				if PLASMID_REGEXP.search(header):
					plasmids += 1
					out.write('>kraken:taxid|45202|OriginalID|{}\n{}\n'.format(header[13:], seq))
				else:
					non_plasmids += 1
					out.write('>{}\n{}\n'.format(header, seq))
		print(plasmids, ' total plamid sequences')
		print(non_plasmids, ' non plasmid sequences')
	
		
		
if __name__ == '__main__':
	kraken_sequences = sys.argv[1]
	output_sequences = sys.argv[2]
	rename_plasmid_seqs(kraken_sequences, output_sequences)
	
#9257  total plamid sequences
#22060  non plasmid sequences

# Run dust on edited_library.fasta
# 

#sed '/^>/! s/[agct]/N/g'  dust_library.fasta	> hard_mask_dust_library.fasta

#kraken-build --add-to-library meglab-kraken-v1-custom-db/edited_vector.fasta --db Custom_kraken_v1_FEB2018/

#kraken-build --add-to-library Custom_kraken_v1_FEB2018/hard_mask_dust_library.fasta --db Custom_kraken_v1_FEB2018/


#kraken-build --build --db Custom_kraken_v1_FEB2018/ --threads 30 --jellyfish-hash-size 7000M
