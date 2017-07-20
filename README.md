# meglab-kraken-custom-db
# Created by MEGlab https://meg.colostate.edu/ 
# Authors: 
  EnriqueDoster and Steven Lakin (lakinsm)
# Description: 
  Python script that downloads and edit bacterial genomes from refseq (reference and representative) to create a custom Kraken database with genomic sequences seperated from plasmid sequences for improved read classification. 
# Usage :
  First must edit the UPPERCASE variables in the kraken-parser.py to make compatible with your machine. 
  Run the script in the directory where you want the script to create a directory with the downloaded genomes:
        
        ./kraken-parser.py new

  Options:
        
        new | Will download the latest "assembly_summary.txt" file from NCBI before downloading genomes


## Workflow
  1. If "new" is used, the script begins by downloading the latest "assembly_summary.txt" file from NCBI before downloading genomes
  2. Else, without the "new" option the script begins by searching the "assembly_summary.txt" for all genomes classified as "representative genome" or "reference genome" and identifying the FTP address, accession number, and taxaid for each.
  3. Then, for each of the FTP addresses the genomic file is downloaded and parsed so that the sequences with headers labeled as plasmids are edited using the string ">kraken:taxid|32630" which kraken will classify as taxon 32630 ("synthetic construct"). All other genomic sequences will be similarly edited ">kraken:taxid|XXXX" to include their corresponding taxaid in their headers.
  4. Using software ,DUST (https://www.ncbi.nlm.nih.gov/pubmed/16796549), each genome is analyzed and low complexity regions of all genomes are identified and nucleotides substituted to lower case nucleotides. (for example: AAATTT => aaattt)
  5. In order to mask the low complexity regions when used with Kraken, the script then parses each sequence to subsitute the lowercase nucleotides to the letter 'N'.
  6. As the script does this for each of the nearly 5,000 genomes, expect the script to take a few hours as the rate limiting step is downloading the genomes from NCBI. 
  
## Once the script is complete, follow the instructions in the Kraken manual to create a custom database   
http://ccb.jhu.edu/software/kraken/MANUAL.html

Briefly, you need to add each of the genomic fasta files using the "kraken-build --add-to-library $file --db $DBNAME" function like shown in the manual:

        for file in output_directory/*.fasta
        do
            kraken-build --add-to-library $file --db $DBNAME
        done

Also, download the taxonomy information:
        
        kraken-build --download-taxonomy --db $DBNAME

Finally, build the database

        kraken-build --build --db $DBNAME

