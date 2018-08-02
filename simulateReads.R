#/usr/bin/Rscript

# Simulates RNA-seq reads using polyester package in R
# fold_changes was edited here to encode DE for specific tRNAs for testing purposes - this can be removed and fold_changes edited to whatever necessary
# fasta should be set to fasta file of transcript sequences as required - here misincorporation mimicking data was generated for testing

library(polyester)
library(Biostrings)

fasta = readDNAStringSet("hg19_misincorporatedSimulatedTranscr.fa")
fold_changes = matrix(rep(1,length(fasta)), nrow = length(fasta))
# encode fold change of 4 (log2FC = 2) for a Gln singleton tRNA (with closely related isodecoders that cluster into a large group)
fold_changes[which(names(fasta) == "Homo_sapiens_tRNA-Gln-CTG-9-1 <unknown description>")] = 4
# encode fold change of 4 for Asn-GTT isoacceptor family of 2 - many other families for same isodecoder present
fold_changes[which(names(fasta) == "Homo_sapiens_tRNA-Asn-GTT-15-1 <unknown description>")] = 4
fold_changes[which(names(fasta) == "Homo_sapiens_tRNA-Asn-GTT-19-1 <unknown description>")] = 4
simulate_experiment("hg19_misincorporatedSimulatedTranscr.fa", readlen = 30, fraglen = 30, fragsd = 0, reads_per_transcript = 1500, paired = FALSE, fold_changes = fold_changes, outdir = "./", num_reps = c(3,3))
