import os
import sys
import re

#req_sam = []
#with open("/home/rewatit/rem_methylation.txt", "r") as samples:
#    for s in samples:
#        s = s.rstrip()
#        req_sam.append(s)

F = os.listdir("/projects/users/rewatit/rna_seq_analysis")
print(F)
for r in F:
    if os.path.isdir("/projects/users/rewatit/rna_seq_analysis/" + r):
        os.chdir("/projects/users/rewatit/rna_seq_analysis/" + r)
        allfiles = os.listdir(".")
        for a in allfiles:
            if a.endswith(".sorted.bam"):
                val = a.split("_")[0]
                os.system("coverageBed -split -abam " + a + " -b /home/rewatit/Rattus_norvegicus_exonic_parts.gff | awk 'BEGIN{OFS=\"\t\"} {print $1,$4,$5,$5-$4+1,$9,$10}' | sort -k 5 > " + val + "_exonic_parts.inclusion")
                #os.system("sed 's/,/\t/g' " + a + " | grep -v description | awk '{OFS=\"\t\"}{print $1,$2+$13, $3-$14,$4,$5,$6}' > " + val + "_intron.bed")
                #os.system("intersectBed -wao -f 1.0 -s -a /home/rewatit/Rattus_norvegicus_exonic_parts.gff -b " + a  + " | awk 'BEGIN{OFS=\"\t\"}{$16 == 0? s[$9] += 0:s[$9] += $14}END{for (i in s) {print i,s[i]}}' | sort -k 1 > " + val  + "_exonic_parts.exclusion")
                #exc = a.split(".")[0] + ".exclusion"
                #os.system("readLength=100")
                #os.system("paste " + a + " " + exc + " | awk -v \"len=100\" 'BEGIN{OFS=\"\t\"; print \"exon_ID\" , \"length\" , \"inclusion\" , \"exclusion\" , \"PSI\"}{NIR=$6/($4+len-1) ; NER=$8/(len-1)}{print $5,$4,$6,$8,(NIR+NER<=0)? \"NA\":NIR / (NIR + NER)}' > exonic_parts.psi")

    os.chdir("/projects/users/rewatit/rna_seq_analysis")
