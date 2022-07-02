#! /usr/bin/env python3

import argparse # to write command-line interfaces
import glob # to retrieve files/pathnames matching a specified pattern
from Bio import SeqIO # for processing sequential data to be fed into downstream sequence models
import re # regular expression
from Bio.SeqUtils import GC # to obtain the GC content
from collections import defaultdict # save a dicc into another dicc
import pandas as pd # data analysis
import subprocess # to import bash functions
import statistics # to calculate statistics


# parser commands:
parser = argparse.ArgumentParser(description = "This program was developed during the MSc in Bioinformatics's thesis of Ester-María López García with the scope of analysing viral communities \
                                                from the metagenomes extracted during the Malaspina Circumnavigation Expedition 2010. \
                                                The pipeline calls many programs in order to perform the assembly of raw reads until the scaffold level, computes statistics about sequences, \
                                                identifies viral scaffolds and discards prokaryotes, executes a QC of scaffolds, calculates relative abundances and standardises to FPKM, \
                                                predicts genes (included Auxiliary Metabolic Genes), conducts taxonomic classification and host prediction and, finaly, annotates genomes. \n\n")

#parser.add_argument("--assemblies", help="File of assemblies", nargs="+", type=str) # file name
parser.add_argument("--assemble", help="Flag to run the assembly using SPAdes for metagenomes. Arguments --r1_files and --r2_files are required.", type=bool, default = False, metavar='')
parser.add_argument("--directory", help="Select a directory from which taking files", type=str, metavar='') # write the directory path
parser.add_argument("--r1_files", help="List of files with forward paired-end reads", type=str, metavar='')
parser.add_argument("--r2_files", help="List of files with reverse paired-end reads", type=str, metavar='')

parser.add_argument("--output_name", help="Write the desired output name", type=str, default=False, metavar='') # pandas output
parser.add_argument("--length", help="Specify the minimum length of the sequences", type=int, default = 999, metavar='')
parser.add_argument("--threads", help="Specify how many threads are needed", type=int, default = 48, metavar='')
parser.add_argument("--filter", help="Filter sequences according to a specified length", type=bool, default = False, metavar='')
parser.add_argument("--concatenate_fasta", help="Run the cat command only when specified", type=bool, default = False, metavar='')
parser.add_argument("--concatenate_matrix", help="Run the cat command only when specified", type=bool, default = False, metavar='')
parser.add_argument("--info_scaffolds", help="Obtain scaffolds information only when specified", type=bool, default = False, metavar='')

parser.add_argument("--run_vibrant", help="Run VIBRANT only when specified", type=bool, default = False, metavar='')
parser.add_argument("--vibrant_output_analysis", help="Run this command only when specified", type=bool, default = False, metavar='')
parser.add_argument("--fasta_merged_vibrant", help="Write the desired VIBRANT merged fasta output file name", type=str, default = "combined_fna_vibrant.fasta", metavar='')
parser.add_argument("--run_prodigal", help="Run Prodigal only when specified", type=bool, default = False, metavar='')
parser.add_argument("--prodigal_output_faa", help="Write the desired Prodigal output FAA file name", type=str, default = "cds.faa", metavar='')
parser.add_argument("--prodigal_output_fna", help="Write the desired Prodigal output FNA file name", type=str, default = "cds.fna", metavar='')
parser.add_argument("--prodigal_output_gff", help="Write the desired Prodigal output GFF file name", type=str, default = "coordinates.gff", metavar='')

# output names
parser.add_argument("--tsv_vibrant", help="Write the desired VIBRANT merged tsv output file name", type=str, default = "table_virus_vibrant.tsv", metavar='') # write the VIBRANT output tsv name
parser.add_argument("--tsv_rafah", help="Write the desired RaFAH merged tsv output file name", type=str, default = "table_virus_rapha.tsv", metavar='')
parser.add_argument("--tsv_checkv", help="Write the desired CheckV merged tsv output file name", type=str, default = "table_virus_checkv.tsv", metavar='')
parser.add_argument("--tsv_vpfclass", help="Write the desired VPFClass merged tsv output file name", type=str, default = "table_virus_vpfclass.tsv", metavar='')
# Classification and host prediction
parser.add_argument("--run_checkV", help="Run checkV only when specified", type=bool, default = False, metavar='')
parser.add_argument("--run_rafah", help="Run Rafah only when specified", type=bool, default = False, metavar='')
parser.add_argument("--run_VPFclass", help="Run VPFclass only when specified", type=bool, default = False, metavar='')
parser.add_argument("--run_bowtie2", help="Run Bowtie2 only when specified", type=bool, default = False, metavar='')
parser.add_argument("--normalise", help="Normalise matrix of counts only when specified", type=bool, default = False, metavar='')
# Annotation
parser.add_argument("--annotation", help="Perform functional annotation with UniProt, KEGG and Pfam databases", type=bool, default = False, metavar='')
parser.add_argument("--uniprot_output", help="Perform functional annotation with UniProt", type=str, default = "output.blastp", metavar='')
parser.add_argument("--kegg_output_hmmsearch", help="Perform functional annotation with KEGG (KOfam)", type=str, default = "kegg_output.hmmsearch", metavar='')
parser.add_argument("--pfam_output_hmmsearch", help="Perform functional annotation with Pfam", type=str, default = "pfam_output.hmmsearch", metavar='')
parser.add_argument("--find_amgs", help="Perform a recovery of the functional annotation information retrieved from Diamond and Hmmsearch (--annotation)", type=bool, default = False, metavar='')

args = parser.parse_args()
parser.print_help()


def main():

    # ------------------------------------- Assembly  ---------------------------------------
    if args.assemble:
        reads_direct = glob.glob(f'{args.directory}*.gz')

        # select reads only from vertical profiles
        samples_code = {} # key: sample name; values: r1 and r2
        samples_code_list = []
        with open("/mnt/lustre/repos/bio/projects/malaspina/release1/00_vertical-profiles-0.22-size-fraction-sample-codes.txt") as file:
            for line in file:
                line = line.strip("\n")
                samples_code[line] = []
                samples_code_list.append(line)

        for elem in reads_direct:
            # file name
            obj_mp_number = re.search("MP(\d)+", elem)
            sample_name = obj_mp_number.group()

            # forward
            if "pair1" in elem and sample_name in samples_code_list:
                samples_code[sample_name].append(elem)
            # reverse
            if "pair2" in elem and sample_name in samples_code_list:
                samples_code[sample_name].append(elem)

        for key in samples_code:
            r1_file = ""
            r2_file = ""
            for value in samples_code[key]:
                if "pair1" in value:
                    r1_file = value
                    print(r1_file)
                if "pair2" in value:
                    r2_file = value

            # call the function
            assembly(threads = args.threads, r1 = r1_file, r2 = r2_file, sample_name = key)


    # ------------------------------- Assembly output selection -----------------------------
    if args.filter:
        # collect all "scaffolds.fasta" files
        sample_direct = glob.glob(f'{args.directory}/Assembly_Group_*')
        for elem in sample_direct:
            sample_files = glob.glob(f'{elem}/scaffolds.fasta')
            ##print(sample_files)
            # select the fasta file
            for file in sample_files:
                # call the function
                processing_fasta(file)


    if args.output_name:
        # create a table from the dictionary
        data_frame = pd.DataFrame.from_dict(dic_seq_info)
        # write the data frame into a new file
        data_frame.index.name = "new_id"
        data_frame.to_csv(args.output_name, sep="\t", na_rep="NA")


    if args.concatenate_fasta:
        # fusion all filtered fasta files into one
        command = "cat filtered_renamed*fasta > filtered_renamed_all.fasta"
        subprocess.call(command, shell=True)

    # ------------------------------- Assembly output information ---------------------------
    if args.info_scaffolds:
        # call the function to extract information about scaffolds (File, Contigs/Scaffolds, Bases, Max, N50, N90)
        scaffolds_direct = f"{args.directory}/*"
        command = f"perl 1_assembly_stats_table.pl {scaffolds_direct}"
        subprocess.call(command, shell = True)

        # calculate average size and sd of gaps for all scaffolds
        scaffolds_files = glob.glob(f"{args.directory}/*")
        gaps_length(scaffolds_files)

    # -------------------------- Viral scaffolds identification ------------------------------
    if args.run_vibrant:
        # call the function to run VIBRANT
        vibrant_function()

    if args.vibrant_output_analysis:
        # call the function to create useful files after with VIBRANT output information
        vibrant_analyse_output()

    # -------------------------- Quality control of viral scaffolds --------------------------
    if args.run_checkV:
        # call the function to run CheckV
        checkv_function()

    # ---------------------------- Scaffolds relative abundance ------------------------------
    if args.run_bowtie2:
        # call the function to run Bowtie2
        bowtie()

    if args.normalise:
        normalisation()

    # ---------------------------- Scaffolds CDS identification ------------------------------
    if args.run_prodigal:
        # call the function to run Prodigal
        prodigal_function()

    # ----------------------------- Viral taxonomy classification ----------------------------
    if args.run_VPFclass:
        # call the function to run VPFclass
        vpfclass_function()

    # ------------------------------------ Host prediction -----------------------------------
    if args.run_rafah:
        # call the function to run RaFAH
        rafah_function()


    if args.concatenate_matrix:
        command = f"Merge_Matrixes.py {args.tsv_checkv} {args.tsv_vpfclass} {args.tsv_rafah}"
        subprocess.call(command, shell=True)

    # ---------------------------- Perform the functional annotation -------------------------
    if args.annotation:
        functional_annot()

    if args.find_amgs:
        amg_hunter()



# function to perform the assembly
def assembly(threads,r1,r2, sample_name):
    print("Assembling reads...", r1, "and", r2)

    # open a new file
    out_assembly = "Assembled_" + sample_name + ".fasta"

    # call SPAdes
    command = f'spades.py --meta --memory 500 --threads {threads} -1 {r1} -2 {r2} -o 00_Assembly_Output/Assembly_Group_{sample_name}' # --meta (assembly metagenomes), --memory (Gb), --threads, -1 (forward reads), -2 (reverse reads), o (output)
    subprocess.call(command, shell=True)


def gaps_length(scaffolds_files):
    print("File\tNumber_gaps\tAverage\tSD\tMax\tMin")
    total_list = []

    for file in scaffolds_files:
        print(f"Processing {file}...")
        # short name
        file_name = file.split("/")[-1]
        eof = False
        with open(file, "r") as IN:
            list_gaps = []
            # read all lines of the file
            while not eof:
                line = IN.readline()
                if not line:
                    eof = True
                # when it finds a GAP
                obj_gap = re.search("NN+",line)
                if obj_gap:
                    # print(obj_gap)
                    gap = obj_gap.group(0)
                    list_gaps.append(len(gap)) # add length to the list
                    total_list.append(len(gap))  # add length to the list

            # once the list is complete, calculate average gap length, sd, max, min by file
            n_gaps = len(list_gaps)
            avr = round(statistics.mean(list_gaps),2)
            sd = round(statistics.stdev(list_gaps),2)
            max_g = max(list_gaps)
            min_g = min(list_gaps)
            print(f"{file_name}\t{n_gaps}\t{avr}\t{sd}\t{max_g}\t{min_g}")

    n_gaps = len(total_list)
    avr = round(statistics.mean(total_list),2)
    sd = round(statistics.stdev(total_list), 2)
    max_g = max(total_list)
    min_g = min(total_list)
    # print("File\tNumber_gaps\tAverage\tSD\tMax\tMin")
    print(f"All\t{n_gaps}\t{avr}\t{sd}\t{max_g}\t{min_g}")



# function to filter fasta files from each sample
def processing_fasta(fasta_file):
    print("Processing...", fasta_file)
    dic_seq_info = defaultdict(dict)
    
    # file name
    obj_mp_number = re.search('MP(\d)+',fasta_file)
    sample_name = obj_mp_number.group(0) # locates the whole match expression 
        
    # open a new file
    out_sample_mil = "filtered_renamed_" + sample_name + ".fasta"
    with open(out_sample_mil,'w', newline='') as OUT:
   
        for seq in SeqIO.parse(fasta_file,'fasta'): # seq=object
                       
            # obtain the sequence id
            seq_id = seq.id
            
            # obtain the sequence new id
            new_id = seq_id.replace("NODE",sample_name)
            new_id = re.sub("_length_(.)+","",new_id) # eliminate from the second _
            
            # change sequence id in the object
            seq.id = new_id
            
            # obtain the sequence
            length = len(seq.seq)        
            ##print(seq_id, new_id, length)
                        
            # GC %
            perc_GC = round(GC(seq.seq),2)
            ##print(seq_id, new_id, length, round(perc_GC,2))
            
            # 3D diccionary
            dic_seq_info["GC"][new_id] = perc_GC
            dic_seq_info["original_id"][new_id] = seq_id 
            dic_seq_info["length"][new_id] = length
            dic_seq_info["original_file"][new_id] = fasta_file
            
            # create a new file with >1000bp sequences
            if length > args.length:
                SeqIO.write(seq, OUT, "fasta")



def vibrant_function():
    # call VIBRANT
    direct_fasta_filtered = glob.glob(f'{args.directory}/filtered_renamed_*fasta')
    for elem in direct_fasta_filtered:
        command_vibrant = f"VIBRANT_run.py -i {elem} -t {args.threads}" # -i (input), -t (threads)
        subprocess.call(command_vibrant, shell=True)


def vibrant_analyse_output():
    # create a 3D diccionary to store all relevant information about sequences from VIBRANT output
    dic_seq_info = defaultdict(dict)

    # after having generated the VIBRANT results
    # create a fasta file with all .fna sequences (complete nucleotide sequences)
    with open(args.fasta_merged_vibrant,'w', newline='') as OUT:

        # select the fna file of each sample
        vibrant_phages_fna = glob.glob(f'VIBRANT_filtered_renamed_MP*/VIBRANT_phages_filtered_renamed_MP*/*.phages_combined.fna')

        # obtain elements from fna file (fasta file)
        for file in vibrant_phages_fna:
            for seq in SeqIO.parse(file,'fasta'): # seq=object

                # obtain the sequence id
                seq_id = seq.id
                # obtain sequence description
                seq_descr = seq.description

                # look for "fragment" in the description in order to know it during downstream analyses
                match = re.search('fragment_(\d)+',seq_descr)
                if match:
                    # change the ID
                    seq.id = seq_id + "_" + match.group()

                # write the output file
                SeqIO.write(seq, OUT, "fasta")


                # create a table with viral sequences information
                # obtain the sequence
                length = len(seq.seq)
                # GC %
                perc_GC = round(GC(seq.seq),2)
                # 3D diccionary
                dic_seq_info["GC"][seq.id] = perc_GC
                dic_seq_info["length"][seq.id] = length
                dic_seq_info["description"][seq.id] = seq_descr.split(' ',1)[1]
                dic_seq_info["AMG_list"][seq.id] = []
                dic_seq_info["type"][seq.id] = ''
                dic_seq_info["Quality"][seq.id] = ''
                dic_seq_info["original_file"][seq.id] = file


    # add a new variable to the dictionary (AMG information: scaffold + AMG KO)
    vibrant_AMG_individuals = glob.glob(f'VIBRANT_filtered_renamed_MP*/VIBRANT_results_filtered_renamed_MP*/VIBRANT_AMG_individuals_filtered_renamed_MP*')
    # open each tsv file
    for tsvfile in vibrant_AMG_individuals:
        # pandas
        info_data_frame = pd.read_csv(tsvfile,sep="\t",index_col="protein",header=0)
        #print(info_data_frame)
        #print(info_data_frame.describe())

        # iterate each row
        for i,row in info_data_frame.iterrows():
            #print(i) # name of index column
            scaffold = row["scaffold"]
            amg = row["AMG KO"]

            # look for "fragment" in the description in order to make future analysis
            match = re.search('fragment_(\d)+',scaffold)
            if match:
                # change the ID
                scaffold = scaffold.split(" ",1)[0] + "_" + match.group()

            else:
                scaffold = row["scaffold"].split(" ",1)[0]

            # create a list into the dictionary with all AMGs contained in each sequence
            dic_seq_info["AMG_list"][scaffold].append(amg)
            
    
    # add a new variables to the dictionary (type: lytic or lysogenic; Quality)
    vibrant_type_quality = glob.glob(f'VIBRANT_filtered_renamed_MP*/VIBRANT_results_filtered_renamed_MP*/VIBRANT_genome_quality_filtered_renamed_MP*')
    # open each tsv file
    for tsvfile in vibrant_type_quality:
        # pandas
        info_data_frame = pd.read_csv(tsvfile,sep="\t",index_col="scaffold",header=0)
        #print(info_data_frame)
        #print(info_data_frame.describe())
        
        # iterate each row
        for i,row in info_data_frame.iterrows():
            #print(i) # name of index column
            type_lyt_lys = row["type"] 
            quality = row["Quality"]
            
            # look for "fragment" in the description in order to make future analysis
            match = re.search('fragment_(\d)+',i)
            if match:
                # change the ID
                identifier = i.split(" ",1)[0] + "_" + match.group() 
                
            else:
                identifier = i.split(" ",1)[0]
            
            # create two new entries in the dictionary with type and quality of each sequence
            dic_seq_info["type"][identifier] = type_lyt_lys
            dic_seq_info["Quality"][identifier] = quality
            
    
    
    # create a table from the dictionary
    data_frame = pd.DataFrame.from_dict(dic_seq_info)
    # write the data frame into a new file 
    data_frame.index.name = "seq_id"
    ###data_frame.to_csv("malaspina_viral_vibrant.tsv",sep="\t",na_rep="NA")
    data_frame.to_csv(args.tsv_vibrant,sep="\t",na_rep="NA")


def checkv_function():
    command = f'checkv end_to_end /mnt/lustre/scratch/elopez/2_VIBRANT_results/combined_fna_vibrant.fasta checkV_output -t {args.threads}'
    subprocess.call(command, shell=True)

    # create empty dictionary for CheckV information
    dic_seq_info = defaultdict(dict)

    # extract valuable information and add to the table (malaspina_viral_table.tsv)
    tsvfile = "/mnt/lustre/scratch/elopez/3_checkV_output/quality_summary.tsv"
    checkv_info = pd.read_csv(tsvfile, sep="\t", index_col="contig_id", header=0)
    # iterate each row
    for i, row in checkv_info.iterrows():
        # create two new entries in the dictionary with completeness and contamination of each sequence
        dic_seq_info["Completeness"][i] = row["completeness"]
        dic_seq_info["Contamination"][i] = row["contamination"]

    # create a table from the dictionary
    data_frame = pd.DataFrame.from_dict(dic_seq_info)
    # write the data frame into a new file
    data_frame.index.name = "seq_id"  # the rest of columns will have their respective names
    data_frame.to_csv(args.tsv_checkv, sep="\t", na_rep="NA")  # "malaspina_virus_checkv.tsv"




# Bowtie2
def bowtie():
    command = f'bowtie2-build /mnt/lustre/scratch/elopez/3_checkV_output/Profiles_Malaspina_ProandViruses_CheckV_Trimmed_Genomes.fasta db_combined' # *.fasta (input), db_combined (output)
    subprocess.call(command, shell=True)  # we call it first, but when we have the db, then we call Bowtie2 for other tasks (make local alignments)

    # table with pair_1 and pair_2 by sample
    csvtable = "/mnt/lustre/scratch/elopez/5_bowtie_results/malaspina_pairs_list.csv"
    # read table with pandas
    pairs_data_frame = pd.read_csv(csvtable, sep=";", index_col="Sample_name", header=0)  # header 0 (there is a header)

    # iterate each row
    for i, row in pairs_data_frame.iterrows():
        # print(i)=sample_name # name of index column
        pair1_file = row["Pair_1"]
        pair2_file = row["Pair_2"]

        outfile = f"{i}xdb_combined"

        command = f"bowtie2 -x db_combined -q -1 {pair1_file} -2 {pair2_file} -S {outfile}.sam --sensitive-local --no-unal --threads {args.threads}"  # db_combined (database prefix), q (fastq), sensitive (local or global), no (only aligned reads), S (SAM)
        subprocess.call(command, shell=True)

        command = f'samtools view -bS {outfile}.sam > {outfile}.bam' # convert SAM file into BAM file
        subprocess.call(command, shell=True)
        command = f'samtools sort {outfile}.bam -o {outfile}.sorted.bam' # sort the BAM file to be able to generate an index
        subprocess.call(command, shell=True)
        command = f'rm -f {outfile}.sam {outfile}.bam' # remove SAM file
        subprocess.call(command, shell=True)
        command = f'samtools index {outfile}.sorted.bam' # create an index from BAM file
        subprocess.call(command, shell=True)
        command = f'samtools idxstats {outfile}.sorted.bam > {outfile}.Counts.tsv' # create a TSV file with the scaffolds name, length and the number of times that the read was mapped over the scaffold
        subprocess.call(command, shell=True)


# Normalise matrix of counts to FPKM
def normalisation():
    length_info = dict() # key: scaffold name; values: length
    map_info = dict() # key: sample; values: total number of reads mapped
    raw_abund_info = defaultdict(dict) # key: sample; values: (key: scaffold, values: number of reads by scaffold)
    fpkm_abund_info = defaultdict(dict) # key: sample; values: (key: scaffold, values: fpkm)

    bowtie_counts_files = glob.glob('*Counts.tsv')

    # count matrix
    for file in bowtie_counts_files:
        print(f"Processing {file}")
        sample = file.split("x")[0]
        with open(file, 'r') as IN:
            for line in IN.readlines():
                vals = line.split('\t')
                if (vals[0] != "*"):
                    (scaffold, length, read_count, dummy) = vals
                    length_info[scaffold] = int(length)
                    raw_abund_info[sample][scaffold] = int(read_count) # for the matrix of counts
                    map_info[sample] = (int(read_count) + map_info.get(sample, 0))


    # FPKM matrix
    for sample in map_info:
        print(f"Calculating FPKM for {sample}")
        for scaffold in length_info:
            # print("Scaffold",scaffold,"Raw abundance",raw_abund_info[sample][scaffold],"Total mapped reads",map_info[sample])
            fpkm_abund_info[sample][scaffold] = ((raw_abund_info[sample][scaffold] * 1000000000) / (map_info[sample] * length_info[scaffold]))

    fpkm_matrix = pd.DataFrame.from_dict(fpkm_abund_info)
    fpkm_matrix.index.name = 'Sequence'
    fpkm_matrix = fpkm_matrix.to_csv("FPKM.tsv", sep="\t", na_rep='NA')



# Prodigal (proteins)
def prodigal_function():
    command = f'prodigal -q -p meta -a {args.prodigal_output_faa} -d {args.prodigal_output_fna} -f gff -i {args.fasta_merged_vibrant} -o {args.prodigal_output_gff}'  # a (aa), d (DNA), f (output format), i (input file), o (output), q (quietly, no stderr), 
    subprocess.call(command, shell=True)



# Taxonomical classification
def vpfclass_function():
    command = f'vpf-class --data-index /mnt/lustre/scratch/elopez/4_VPFClass_output/data-index.yaml --input-seqs /mnt/lustre/scratch/elopez/4_checkV_output/Prodigal_output_2/Profiles_Malaspina_ProandViruses_CheckV_Trimmed_Genomes.fasta -o /mnt/lustre/scratch/elopez/4_VPFClass_output/After_checkV_output --chunk-size 1000 --workers {args.threads}'
    subprocess.call(command, shell=True)

    # create empty dictionary for VPFClass output information
    dic_seq_info = defaultdict(dict)

    # extract valuable information and add to the table (malaspina_viral_table.tsv)
    # select the fna file of each sample
    vpfclass_folder = glob.glob(f'./4_VPFClass_output/After_checkV_output/*tsv')

    # obtain elements from fna file (fasta file)
    for file in vpfclass_folder:
        # print(file)
        # files: baltimore.tsv, family.tsv and genus.tsv
        match = re.search(f'host*', file)
        if not match:
            vpfclass_info = pd.read_csv(file, sep="\t", index_col="virus_name", header=0)

            # search_string_start = re.compile("./4_VPFClass_output/")
            search_string_start = re.compile("./4_VPFClass_output/After_checkV_output/")
            search_string = re.compile("\.tsv")  # interprets only one time

            class_name = re.sub(search_string_start, "", file) # family.tsv
            class_name = re.sub(search_string, "", class_name) # family
            memb_ratio = re.sub(search_string_start, "", file) # family.tsv
            memb_ratio = re.sub(search_string, "_membership", memb_ratio) # family_membership

            print(f'Processing... {file}')
            # iterate each row
            for i, row in vpfclass_info.iterrows():

                # get the higher membership_ratio
                membership_value = dic_seq_info[memb_ratio].get(i, 0)  # if there is no value, gives 0. Otherwhise, it gives the value
                # print(membership_value, row["membership_ratio"])

                print(i, row)
                if row["membership_ratio"] > membership_value:
                    # create one new entry in the dictionary with class_name of each sequence
                    dic_seq_info[class_name][i] = row["class_name"]
                    dic_seq_info[memb_ratio][i] = row["membership_ratio"]

            # remove membership_ratio key from the dictionary
            del dic_seq_info[memb_ratio]

    # create a table from the dictionary
    data_frame = pd.DataFrame.from_dict(dic_seq_info)
    # write the data frame into a new file
    data_frame.index.name = "seq_id"  # the rest of columns will have their respective names
    data_frame.to_csv(args.tsv_vpfclass, sep="\t", na_rep="NA")  # "malaspina_virus_checkv.tsv" / 2nd time: "malaspina_trimmed_checkV_virus_vpfclass.tsv"




# Host prediction
def rafah_function():
    command = f'./4_RaFAH_output/After_checkV_output/RaFAH_v0.2.pl --predict --merged_cds_file_name /mnt/lustre/scratch/elopez/4_checkV_output/Prodigal_output_2/Profiles_Malaspina_ProandViruses_CheckV_Trimmed_Genomes.faa --min_cutoff 0 --threads {args.threads} --file_prefix RaFAH_malaspina_'
    subprocess.call(command, shell=True)

    # create empty dictionary for RaFAH output information
    dic_seq_info = defaultdict(dict)
    
    # extract valuable information and add to the table (malaspina_viral_table.tsv)
    #tsvfile = "/mnt/lustre/scratch/elopez/4_RaFAH_output/RaFAH_malaspina__Seq_Info_Prediction.tsv"
    tsvfile = "/mnt/lustre/scratch/elopez/4_RaFAH_output/After_checkV_output/RaFAH_malaspina__Seq_Info_Prediction.tsv"
    rafah_info = pd.read_csv(tsvfile,sep="\t",index_col="Variable",header=0)
    # iterate each row
    for i,row in rafah_info.iterrows():        
        # create two new entries in the dictionary with the predicted host and the predicted host score of each sequence
        dic_seq_info["Predicted_Host"][i] = row["Predicted_Host"] 
        dic_seq_info["Predicted_Host_Score"][i] = row["Predicted_Host_Score"]
    
    # create a table from the dictionary
    data_frame = pd.DataFrame.from_dict(dic_seq_info)
    # write the data frame into a new file 
    data_frame.index.name = "seq_id" # the rest of columns will have their respective names 
    data_frame.to_csv(args.tsv_rafah,sep="\t",na_rep="NA") # "malaspina_virus_rapha.tsv" / 2nd time: "malaspina_trimmed_checkV_virus_rapha.tsv"



# Functional annotation
def functional_annot():
    
    # UniProt 
    db_file = "/mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.dmnd"
    command = f'diamond blastp --matrix BLOSUM45 --threads {args.threads} --more-sensitive --db {db_file} --outfmt 6 --query {args.prodigal_output_faa} --max-target-seqs 100 --evalue 0.001 --out {args.uniprot_output}' # db (Uniprot), query (prot sequences to compare. faa), max.target.seqs (top100)
    subprocess.call(command, shell=True)    # Blosum45 distant homology (our case)
    
    # KEGG
    db_file = "/mnt/lustre/scratch/elopez/KOfam/All_KOs.hmm"
    command = f'hmmsearch -o {args.kegg_output_hmmsearch} --noali --cpu {args.threads} {db_file} {args.prodigal_output_faa}'
    subprocess.call(command, shell=True)
    
    # Pfam 
    db_file = "/mnt/lustre/repos/bio/databases/public/pfam/pfam_release_34.0/Pfam-A.hmm"
    command = f'hmmsearch -o {args.pfam_output_hmmsearch} --noali --cpu {args.threads} {db_file} {args.prodigal_output_faa}' # hmmscan and hmmsearch (faster)
    subprocess.call(command, shell=True)


def amg_hunter():
    db_file = "/mnt/lustre/scratch/elopez/KOfam/All_KOs.hmm"
    command = f'python3 AMG_Hunter.py --parse_only True --abundance /mnt/lustre/scratch/elopez/5_bowtie_results/After_checkV_output/FPKM.tsv --cds {args.prodigal_output_faa} --kegg_db {db_file} --kegg_json /mnt/lustre/scratch/elopez/KOfam/ko00001.json --ko_hmm_dir /mnt/lustre/scratch/elopez/KOfam/profiles/'
    subprocess.call(command, shell=True)


main()
