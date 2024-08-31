# For classifying the 16s rRNA gene sequences to taxonomy
# Need to download and add the path of these databases
silva_template=Model_data/Silva_ref_data/silva.nr_v138_1.align
taxo_ref=Model_data/Silva_ref_data/silva.nr_v138_1.tax

# classify every OTU to taxonomy
otu_seqs=$1
mothur "#classify.seqs(fasta=$otu_seqs,reference=$silva_template,taxonomy=$taxo_ref,cutoff=50,processors=16)"