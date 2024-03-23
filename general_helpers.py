import glob, gzip, math, random
from Bio import SeqIO, Align, Entrez
import os
import gzip
import shutil

Entrez.email = "ghassan.abboud@epfl.ch"

def extract_gz(src_file, dst_file):
    """
    Extracts a gzipped file.

    Args:
        src_file (str): Path to the source gzipped file.
        dst_file (str): Path to the destination file.

    Returns:
        None
    """
    with gzip.open(src_file, 'rb') as f_in:
        with open(dst_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def read_fastq(fastq_filepath=None):
    """
    Read a fastq file and return the reads.
    :param fastq_filepath: filepath of the .fastq file [str]
    :return: reads from the fastq file [list of Bio.SeqRecord.SeqRecord]
    """
    if fastq_filepath is None: fastq_filepath = "data/rbcL_Qiagen_tomato.fastq" # default path (example)
    if fastq_filepath.lower().endswith('.gz'):
        f = gzip.open(fastq_filepath, 'rt')
    else:
        f = open(fastq_filepath, 'rt')
    reads = []
    for read in SeqIO.parse(f, "fastq"):
        reads.append(read)
    return reads

def concatenate_fastq(src_folder=None, dst=None):
    """
    Concatenate all .fastq from a folder into a single .fastq file.
    :param folder: folder containing the .fastq files (usually fastq_pass folder) [str]
    :param dst: destination file (.fastq) [str]
    :return: None
    """
    if src_folder is None: src_folder = "fastq_pass" # default folder (example)
    if dst is None: dst = f"{src_folder}/concatenation.fastq" # default destination (example)
    if dst[-6:] != ".fastq": dst = f"{dst}.fastq"
    def get_file_iterator(filename):
        if filename.lower().endswith('.gz'):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'rt')
        return f
    fastq_list = [fastq for fastq in glob.glob(f"{src_folder}/*.fastq*")]
    fastq_iterators = [SeqIO.parse(get_file_iterator(fastq), "fastq") for fastq in fastq_list]
    while True:
        for fq in fastq_iterators:
            try:
                SeqIO.write(next(fq), open(dst,"at"), "fastq")
            except StopIteration:
                fastq_iterators.remove(fq)
        if len(fastq_iterators) == 0:
            break

def download_sequence(species, gene_name, dst, start_length=None, stop_length= None, id = None, permissive_search = True):
    """
    download sequence from GenBank through the Entrez database. 

    Parameters:
    ----------
    species(str): name of species
    gene_name(str): name of gene
    dst(str,Path-like): destination file path
    start_length(int): minimum length of sequence
    stop_length(int): maximum length of sequence
    id(list): list of NCBi ids of sequences to download. If provided, overrides gene_name and species.
    permissive_search(bool, default = True): when True, if Advanced NCBI query returns nothing, replace it with a less precise general query.
    """
    
    if id == None:
        search_term = f"{gene_name}[Gene Name] AND {species}[Organism]"
        if start_length!= None or stop_length!= None:
                search_term += f" {start_length}:{stop_length}[Sequence Length]"
        handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
        search_result = Entrez.read(handle)
        handle.close()
        id = search_result["IdList"]
        n=0
        for i in id:
            n+=1
        if n==0 and permissive_search:
            search_term = f"{gene_name} {species} {start_length}:{stop_length}[Sequence Length]"
            handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
            search_result = Entrez.read(handle)
            handle.close()
            id = search_result["IdList"]
    n=0
    for i in id:
        n+=1
    if n==1:
        handle = Entrez.efetch(db="nucleotide", id=id, retmode = "fasta", rettype = "fasta")
        sequence = handle.read()
        handle.close()
        with open(dst, mode="a") as writer:
            writer.write(sequence)
    
def extract_fastq(main_dir, sample_nb, dst):
    """
    extract all fastq files under a certain barcode/sample number from an expedition results folder.

    Parameters
    ----------
    main_dir: expedition folder path
    sample_nb: sample number for which to extract fastq. ie 6 to extract from folder barcode06
    dst: destination file path
    """
    for root, dirs, files in os.walk(main_dir):
        if "fastq_pass" in dirs:
            pass_dir = ospath.join(root, "fastq_pass")
            break
    for root, dirs, files in os.walk(pass_dir):
        if f"barcode{sample_nb}" in dirs:
            fastq_dir = ospath.join(root, f"barcode{sample_nb}")
        elif f"barcode0{sample_nb}" in dirs:
            fastq_dir = ospath.join(root, f"barcode0{sample_nb}")
    concatenate_fastq(fastq_dir, dst)