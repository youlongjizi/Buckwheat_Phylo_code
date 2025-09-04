import re
import os
import argparse

def process_fasta_to_phy(input_file, output_file):
    """
    Converts a FASTA file to PHYLIP format.

    Args:
    input_file (str): Path to the input FASTA file.
    output_file (str): Path to the output PHYLIP file.
    """
    # Read and parse the FASTA file
    try:
        with open(input_file, 'r') as fin:
            content = fin.read()
            sequences = re.findall(r'(?m)^>([^ \n]+)[^\n]*\n([^>]*)', content)
    except FileNotFoundError:
        print(f"Error: The file {input_file} does not exist.")
        return
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        return

    # Write to PHYLIP format
    try:
        with open(output_file, 'w') as fout:
            num_sequences = len(sequences)
            if num_sequences > 0:
                sequence_length = len(sequences[0][1].replace("\n", ""))  # Remove newlines from sequence
                fout.write(f'{num_sequences} {sequence_length}\n')
                for identifier, sequence in sequences:
                    # Remove newlines and write sequence in a single line
                    sequence_data = sequence.replace("\n", "")
                    fout.write(f'{identifier}    {sequence_data}\n')
            else:
                print(f"No sequences found in {input_file}")
    except Exception as e:
        print(f"Error writing to {output_file}: {e}")

def batch_process_files(input_dir, output_dir, file_suffix):
    """
    Batch processes files in the input directory and converts them to PHYLIP format.

    Args:
    input_dir (str): Directory containing the input files.
    output_dir (str): Directory to save the output PHYLIP files.
    file_suffix (str): Suffix of the files to process.
    """
    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Process each file with the specified suffix in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(file_suffix):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename.replace(file_suffix, '.phy'))
            process_fasta_to_phy(input_file, output_file)
            if os.path.exists(output_file):
                print(f'Processed {input_file} -> {output_file}')

def main():
    parser = argparse.ArgumentParser(description='Convert FASTA files to PHYLIP format.')
    parser.add_argument('-i', '--input_file', type=str, help='Single input FASTA file')
    parser.add_argument('-o', '--output_file', type=str, help='Single output PHYLIP file')
    parser.add_argument('-d', '--input_directory', type=str, help='Directory containing the FASTA files')
    parser.add_argument('-O', '--output_directory', type=str, help='Directory to save the .phy files')
    parser.add_argument('-s', '--suffix', type=str, default='aln', help='File suffix to process for batch mode (default: aln)')

    args = parser.parse_args()

    if args.input_file and args.output_file:
        # Process a single file
        process_fasta_to_phy(args.input_file, args.output_file)
    elif args.input_directory and args.output_directory:
        # Batch process files in a directory
        batch_process_files(args.input_directory, args.output_directory, '.' + args.suffix)
    else:
        parser.error("Either specify -i and -o for single file processing or -d and -O for batch processing.")

if __name__ == '__main__':
    main()

