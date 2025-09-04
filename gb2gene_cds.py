import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import warnings
from Bio.GenBank import BiopythonParserWarning

# 忽略BiopythonParserWarning警告
warnings.filterwarnings("ignore", category=BiopythonParserWarning)

def extract_sequences(gb_path, output_file, gene_name=None, sequence_type='CDS'):
    try:
        gb_records = SeqIO.parse(gb_path, "genbank")
    except FileNotFoundError:
        raise FileNotFoundError(f"The file {gb_path} does not exist")
    
    extracted_sequences = []
    
    for record in gb_records:
        for feature in record.features:
            if feature.type == sequence_type:
                if gene_name is None or ("gene" in feature.qualifiers and gene_name in feature.qualifiers["gene"]) or ("product" in feature.qualifiers and gene_name in feature.qualifiers["product"]):
                    seq = record.seq[feature.location.start:feature.location.end] if feature.location.strand == 1 else record.seq[feature.location.start:feature.location.end].reverse_complement()
                    seq_record = SeqRecord(seq, id=f"{record.id}_{feature.qualifiers.get('gene', [''])[0]}_{feature.qualifiers.get('product', [''])[0]}", description=f"{record.description} {feature.qualifiers.get('gene', [''])[0]} {feature.qualifiers.get('product', [''])[0]}")
                    extracted_sequences.append(seq_record)
    
    # Write the extracted sequences to a FASTA file
    SeqIO.write(extracted_sequences, output_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Extract genes or CDS sequences from a GenBank file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the GenBank file.")
    parser.add_argument("-o", "--output", required=True, help="Output file name for the extracted sequences.")
    parser.add_argument("-g", "--gene", help="Name of the gene to extract. If not provided, all genes will be extracted.")
    parser.add_argument("-t", "--type", choices=['CDS', 'gene'], default='CDS', help="Type of sequence to extract: CDS or gene.")
    args = parser.parse_args()
    
    extract_sequences(args.input, args.output, gene_name=args.gene, sequence_type=args.type)

if __name__ == "__main__":
    main()

