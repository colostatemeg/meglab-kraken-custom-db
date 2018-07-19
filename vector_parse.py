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
			vectors = 0
			for header, seq in data.items():
				vectors += 1
				out.write('>kraken:taxid|32630|{}\n{}\n'.format(header, seq))
		print(vectors, ' total vector sequences')

	
		
		
if __name__ == '__main__':
	input_sequences = sys.argv[1]
	output_sequences = sys.argv[2]
	rename_plasmid_seqs(input_sequences, output_sequences)
	
# some empty header lines, use the following awk code to only keep useful headers 
#awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' edited_vector.fasta > vector_contaminant.fasta