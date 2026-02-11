
import sys
import os
from collections import defaultdict

def parse_gtf(gtf_file):
    """
    Parses a GTF file and groups exons by transcript_id.
    """
    transcripts = defaultdict(list)
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            
            # We only care about exons for BED12 structure
            if feature_type != 'exon':
                continue
                
            chrom = parts[0]
            start = int(parts[3]) - 1  # 0-based start
            end = int(parts[4])        # 1-based end (exclusive in BED)
            strand = parts[6]
            
            # Handle malformed GTF where attributes might contain tabs
            # specific to this user's viral hybrid GTF
            attributes = " ".join(parts[8:])
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                
                # Handle space or = separation
                if '"' in attr:
                    try:
                        key, val = attr.split('"', 1)
                        val = val[:-1] # Remove trailing quote
                    except ValueError:
                        continue
                else:
                    try:
                        key, val = attr.split(None, 1)
                    except ValueError:
                        continue
                
                attr_dict[key.strip()] = val
            
            transcript_id = attr_dict.get('transcript_id')
            gene_id = attr_dict.get('gene_id')
            
            # Fallback for viral genes that might default to gene_id as transcript
            if not transcript_id:
                transcript_id = gene_id
                
            if not transcript_id:
                continue
                
            transcripts[transcript_id].append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id
            })
            
    return transcripts

def convert_to_bed12(transcripts):
    """
    Converts grouped transcript data to BED12 format.
    """
    for tx_id, exons in transcripts.items():
        if not exons:
            continue
            
        # Sort exons by start position
        exons.sort(key=lambda x: x['start'])
        
        chrom = exons[0]['chrom']
        strand = exons[0]['strand']
        gene_id = exons[0]['gene_id']
        
        # Transcript start and end
        tx_start = exons[0]['start']
        tx_end = exons[-1]['end']
        
        # Block info
        block_count = len(exons)
        block_sizes = []
        block_starts = []
        
        for exon in exons:
            block_sizes.append(str(exon['end'] - exon['start']))
            block_starts.append(str(exon['start'] - tx_start))
            
        # BED12 fields
        # 1. chrom
        # 2. chromStart
        # 3. chromEnd
        # 4. name (transcript_id) 
        # 5. score (0)
        # 6. strand
        # 7. thickStart (CDS start - using tx_start here as simplified)
        # 8. thickEnd (CDS end - using tx_end here as simplified)
        # 9. itemRgb (0)
        # 10. blockCount
        # 11. blockSizes
        # 12. blockStarts
        
        bed_line = [
            chrom,
            str(tx_start),
            str(tx_end),
            tx_id, 
            "0",
            strand,
            str(tx_start),
            str(tx_end),
            "0",
            str(block_count),
            ",".join(block_sizes),
            ",".join(block_starts)
        ]
        
        print("\t".join(bed_line))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python gtf2bed.py <input.gtf>\n")
        sys.exit(1)
        
    gtf_file = sys.argv[1]
    transcripts = parse_gtf(gtf_file)
    convert_to_bed12(transcripts)
