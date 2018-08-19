'''
replace text through files
'''
import os
directory="/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS"
for filename in os.listdir(directory):
    if filename.endswith("allgenes_nexus.fasta"):
        print(os.path.join(directory, filename))
        with open(directory+"/"+filename) as open_file, open(directory+"/SGS_allgenes_final.fasta","w") as new_file:
            linecount=0
            for line1 in open_file:
                linecount+=1 # TEXT TO BE MODIFIED BELOW
                if linecount < 7:
                    continue
                elif ";" in line1:
                    continue
                elif "-" not in line1:
                    continue
                else:
                    line2 = '\n'.join(line1.split())
                    line3 = line2.replace("(","")
                    line4 = line3.replace(")","")
                    if linecount is 7:
                    	new_file.write(">"+line4)
                    else:
                    	new_file.write("\n"+">"+line4)
        continue
    else:
        continue
