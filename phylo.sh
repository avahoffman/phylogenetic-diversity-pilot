#!/bin/bash
# Author:       Ava Hoffman
# Email:        ava.hoffman@colostate.edu
# Date:         2018-5-22
# Description:
# 
#
#

# /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyloGeneratorMac_1-3-rc/phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodSimple -wd /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/examples/ -email avamariehoffman@gmail.com -gene rbcL -species /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyloGeneratorMac_1-3-rc/Demos/Silwood_Plants/species.txt -alignment mafft -phylogen beast-GTR-GAMMA
# 
# #Simple search with dated constraint tree (NOTE: using downloaded sequences now)
# /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyloGeneratorMac_1-3-rc/phyloGenerator.app/Contents/MacOS/phyloGenerator -name silwoodConstrained -wd /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/examples/ -email avamariehoffman@gmail.com -gene rbcL -dna /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyloGeneratorMac_1-3-rc/Demos/Silwood_Plants/sequences.fasta -alignment mafft -phylogen beast-GTR-GAMMA -consTree /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyloGeneratorMac_1-3-rc/Demos/Silwood_Plants/constraint.tre
# 
# # God knows why, but need to install JDK not other Java ?? : http://www.oracle.com/technetwork/java/javase/downloads/index.html
# CLASSPATH=${CLASSPATH}:/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyloGeneratorMac_1-3-rc/phyloGenerator.app/Contents/Resources/requires/beast.jar
# PATH=$PATH:

## MUST DELETE OLD APPLET AND REINSTALL JAVA...
# https://it.uoregon.edu/reverting-java-mac

#
#
#
#

## phylogeny instructions
## make sure .fasta is correctly formatted and species names are the same
## remove any species missing sequences
## align genes using MAFFT : https://www.ebi.ac.uk/Tools/msa/
##		can play around with MAFFT options but not really necessary
##		If there is one alignment that's particularly gap-y, especially in protein coding gene, check if the reverse complement is better:
##		http://www.cellbiol.com/scripts/complement/dna_sequence_reverse_complement.php
## save output as a new .fasta
## Use MEGA https://www.megasoftware.net/ to visualize your aligned sequences
## 		check for gaps of 3 in coding sequences, also use translate portion to make sure things make sense
## Use GBlocks to trim extra stuff off ends: http://molevol.cmima.csic.es/castresana/Gblocks_server.html
## Fill in any missing species with length of trimmed gene.. use dashes (easy to do in excel with =rept("-",701) for example)
## Use phyutility : https://code.google.com/archive/p/phyutility/

## java -jar /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyutility/phyutility.jar -concat -in /Users/avahoffman/Desktop/plot_level_fasta_sgs/A1_ITS_aligned.fasta /Users/avahoffman/Desktop/plot_level_fasta_sgs/A1_matK_aligned_trimmed.fasta /Users/avahoffman/Desktop/plot_level_fasta_sgs/A1_ndhF_aligned_trimmed.fasta /Users/avahoffman/Desktop/plot_level_fasta_sgs/A1_rbcL_aligned_trimmed.fasta -out A1_allgenes.fasta

## java -jar /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyutility/phyutility.jar -concat -in /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS/COMMUNITY_ITS_sequences_SGS_aligned.fasta /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS/COMMUNITY_rbcL_sequences_SGS_aligned_trimmed.fasta /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS/COMMUNITY_TrnL_sequences_SGS_aligned.fasta /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS/COMMUNITY_matK_sequences_SGS_aligned_trimmed.fasta /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS/COMMUNITY_TrnL-TrnFspacer_sequences_SGS_aligned_trimmed.fasta /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS/COMMUNITY_psbA_sequences_SGS_aligned.fasta /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS/COMMUNITY_ndhF_sequences_SGS_aligned_trimmed.fasta -out SGS_allgenes_nexus.fasta

## java -jar /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyutility/phyutility.jar -concat -in ITS_sequences_KNZ_aligned.fasta TrnL-TrnFspacer_sequences_KNZ_aligned.fasta TrnL_sequences_KNZ_aligned_trimmed.fasta matK_sequences_KNZ_aligned_trimmed.fasta ndhF_sequences_KNZ_aligned_trimmed.fasta psbA_sequences_KNZ_aligned_trimmed.fasta rbcL_sequences_KNZ_aligned_trimmed.fasta -out KNZ_allgenes_nexus.fasta

## java -jar /Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/phyutility/phyutility.jar -concat -in

## outputs a nexus file, may need to reformat for other applications.

## may use CIPRES and RAxML to develop best final tree..