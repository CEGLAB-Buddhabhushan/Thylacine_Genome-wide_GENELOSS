import pandas as pd
import sys

# Check if transcript_id is provided as an argument
if len(sys.argv) != 2:
    print("Usage: python codon_pos_bed6.py <transcript_id>")
    sys.exit(1)

# Get the transcript_id from the command-line argument
transcript_id = sys.argv[1]

# Load the GTF file
gtf_file = '/media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/GCF_902635505.1_mSarHar1.11_genomic.gtf.cds.gtf'
gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

# Assign column names
gtf_df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# Extract transcript_id from the attribute column
gtf_df['transcript_id'] = gtf_df['attribute'].str.extract('transcript_id "([^"]+)"')

# Filter the DataFrame for the specific transcript_id
transcript_df = gtf_df[(gtf_df['transcript_id'] == transcript_id) & (gtf_df['feature'] == 'CDS')]

# Initialize a list to store BED6 format codon positions
bed6_entries = []

# Initialize a variable to track the codon phase shift
current_phase_shift = 0

# Initialize codon counter
codon_counter = 1

# Iterate over each CDS segment and calculate codon positions considering the phase
for _, row in transcript_df.iterrows():
    seqname = row['seqname']
    start = row['start']
    end = row['end']
    strand = row['strand']
    phase = int(row['frame'])
    
    # Adjust start position based on the current phase shift
    if strand == '+':
        start += current_phase_shift
        codon_start_positions = range(start, end + 1, 3)
    else:  # For negative strand
        end -= current_phase_shift
        codon_start_positions = range(end, start - 1, -3)
    
    # Update phase shift for the next CDS
    current_phase_shift = (3 - (end - start + 1) % 3) % 3
    
    # Create BED6 entries
    for codon_start in codon_start_positions:
        if strand == '+':
            chromStart = codon_start - 1  # BED is 0-based
            chromEnd = codon_start + 3    # codon is 3 bases long
        else:
            chromStart = codon_start - 3
            chromEnd = codon_start

        # Label the codon sequentially
        codon_name = f"codon{codon_counter}"
        codon_counter += 1
        
        bed6_entries.append([seqname, chromStart, chromEnd, codon_name, 0, strand])

# Convert to DataFrame for easier handling
bed6_df = pd.DataFrame(bed6_entries, columns=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])

# Save the codon positions in BED6 format to a file
output_file = f'{transcript_id}.codon_positions.bed'
bed6_df.to_csv(output_file, sep='\t', header=False, index=False)

print(f"Output saved to {output_file}")
