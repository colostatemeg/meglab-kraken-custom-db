# meglab-kraken-custom-db
# Created by MEGlab https://meg.colostate.edu/ 
# Authors: 
  EnriqueDoster and Steven Lakin (lakinsm)
# Description: 
  Python script that edits the bacterial genomes from refseq (reference and representative) to create a custom Kraken database (version 1) with genomic sequences seperated from plasmid sequences for improved read classification. 

# Workflow :
  First, must follow instructions from Kraken manual to download desired genomes. For example, can download genomes for Archaea, Bacteria, and Viruses.
        
        kraken-build --download-library archaea --db $DBNAME
        
        kraken-build --download-library bacteria --db $DBNAME
        
        kraken-build --download-library viral --db $DBNAME
  Then, concatenate all downloaded genomes in the /library directories. 
        
        cat $DBNAME/library/*/library.fna > $DBNAME/library/concatenated_genomes.fna
  Now, must edit the UPPERCASE variables in the kraken-v1-parser.py to make compatible with your machine. 
  Run the kraken-v1-parser.py script to parse the sequences with headers labeled as plasmids and edited using the string ">kraken:taxid|45202" which kraken will classify as taxon 45202 ("plasmid"). Also, concantenate the "edited_vector.fasta" to your genomes which  includes contaminant vector sequences from: UniVec, EmVec, and the Enterobacteria phage PhiX-174 genome. 
        
        python meglab-kraken-custom-db/kraken-parser/kraken-v1-parser.py $DBNAME/library/concatenated_genomes.fna $DBNAME/library/edited_library.fasta
        
        cat $DBNAME/library/edited_library.fasta meglab-kraken-custom-db/edited_vector.fasta > $DBNAME/library/full_library.fasta
  The resulting file must now be processed with DUST (https://www.ncbi.nlm.nih.gov/pubmed/16796549), each genome is analyzed and low complexity regions of all genomes are identified and nucleotides substituted to lower case nucleotides. (for example: AAATTT => aaattt)

        dustmasker -in $DBNAME/library/full_library.fasta -out $DBNAME/library/dust_full_library.fasta -infmt fasta -outfmt fasta
  In order to mask the low complexity regions when used with Kraken, the lowercase nucleotides need to be substituted to the letter 'N'.

        sed '/^>/! s/[agct]/N/g' $DBNAME/library/dust_full_library.fasta > $DBNAME/mask_dust_full_library.fasta
  
## Once the steps above are complete, follow the instructions in the Kraken manual to build the custom database
http://ccb.jhu.edu/software/kraken/MANUAL.html

Briefly, you need to add the fasta file using the "kraken-build --add-to-library $file --db $DBNAME" function like shown in the manual:

        kraken-build --add-to-library $DBNAME/mask_dust_full_library.fasta --db $DBNAME

Also, download the taxonomy information:
        
        kraken-build --download-taxonomy --db $DBNAME

Finally, build the database

        kraken-build --build --db $DBNAME

