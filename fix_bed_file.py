original_bed_file = open("/Users/m006703/Illumina/Knowledge_Network/files/lifted_over.bed", "r")
fixed_bed_file = open("/Users/m006703/Illumina/Knowledge_Network/files/fixed_hg38_nmpan.bed", "w")

for line in original_bed_file:
    line = line.rstrip()
    line_item = line.split("\t")
    chromosome = line_item[0]
    start = line_item[1]
    stop = line_item[2]
    fixed_start = str(int(start) - 1)
    fixed_stop = str(int(stop) - 1)
    fixed_bed_file.write(chromosome + "\t" + fixed_start + "\t" + fixed_stop + "\n")

original_bed_file.close()
fixed_bed_file.close()
