import os
import re
import argparse
from collections import defaultdict
from collections import namedtuple
from Bio import SeqIO

import pdb

# -gpff
# -augustus_hints
# -best_hits

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    parser = argparse.ArgumentParser(
            description="""Add annotations to a Braker gtf file from a BLAST search of proteins"""
        )
    parser.add_argument(
            '--augustus-hints',
            required=True,
            type=is_file,
            action=FullPaths,
            help='The augustus.hints.gtf file from braker'
        )
    parser.add_argument(
            "--ncbi-gpff",
            required=True,
            type=is_file,
            action=FullPaths,
            help="""The GPFF file from the assembly used to add annotation from Genbank."""
        )
    parser.add_argument(
            "--best-hits",
            required=True,
            type=is_file,
            action=FullPaths,
            help="""The tab-formatted output file from BLAST for the search against proteins."""
        )
    parser.add_argument(
            "--output",
            required=True,
            action=FullPaths,
            help="""The output GTF-formatted file with annotation data."""
        )
    return parser.parse_args()

def main():
    args = get_args()
    print("Parsing the GPFF file")
    genbank = {}
    GbRecord = namedtuple('GbRecord', ['product', 'organism', 'gene', 'db_xref'])
    # parse the genbank file - this just digests information about the
    # annotation in the reference set of proteins we're using to attach
    # names to our transcripts/genes
    with open(args.ncbi_gpff) as infile:
        for record in SeqIO.parse(infile, "genbank"):
            for feat in record.features:
                if feat.type == 'CDS':
                    gene = feat.qualifiers['gene'][0]
                    db_xref = feat.qualifiers['db_xref']
                elif feat.type == 'Protein':
                    product = feat.qualifiers['product'][0]
                elif feat.type == 'source':
                    organism = feat.qualifiers['organism'][0]
            genbank[record.id] = GbRecord(product=product, organism=organism, gene=gene, db_xref=db_xref)
    # this parses the GTF file produced by braker to get a list of the
    # auto-assigned gene and trascript names
    print("Parsing the Hints GTF file")
    gtf_annotation = defaultdict(list)
    with open(args.augustus_hints) as infile:
        for line in infile:
            ls = line.strip().split("\t")
            if ls[2] == 'transcript':
                gene, transcript = ls[-1].split(".")
                gtf_annotation[gene].append(transcript)
    blast_genes = defaultdict(set)
    with open(args.best_hits) as infile:
        for line in infile:
            ls = line.strip().split("\t")
            blast_genes[ls[0]].add(ls[1].split("|")[1])
    # report any gene.transcripts that hit multiple accessions
    for gene, hits in blast_genes.items():
        if len(hits) > 1:
            print("WARN: {} has multiple hits".format(gene))
    print("Parsing the BLAST data")
    annotation = defaultdict(dict)
    annotation_counts = {'annotated':0, 'missing':0}
    for gene, transcripts in gtf_annotation.items():
        for transcript in transcripts:
            accession = blast_genes["{}.{}".format(gene, transcript)]
            if accession != set():
                accession = list(accession)[0]
                annotation[gene][transcript] = genbank[accession]
                annotation_counts['annotated'] += 1
            else:
                annotation[gene][transcript] = None
                annotation_counts['missing'] += 1
    with open(args.output, "w") as outfile:
        with open(args.augustus_hints) as infile:
            for line in infile:
                if line.startswith("#"):
                    outfile.write
                else:
                    ls = line.strip().split("\t")
                    cat, desc = ls[2], ls[-1]
                    if cat == "gene":
                        print(cat, desc)
                        # check to see if this exists in annotation
                        gene = desc
                        if annotation[gene]['t1'] != None:
                            new_desc = 'gene_id "{}"; gene_name "{}"; gene_product "{}"; gene_organism "{}";'.format(
                                desc,
                                annotation[desc]['t1'].gene,
                                annotation[desc]['t1'].product,
                                annotation[desc]['t1'].organism
                            )
                            if annotation[gene]['t1'].db_xref:
                                db_xref = "; ".join(['ncbi_gene_identifier_{} "{}"'.format(i, j) for i,j in enumerate(annotation[gene]['t1'].db_xref)])
                                new_desc = new_desc + db_xref
                        else:
                            new_desc = 'gene_id "{}"'.format(desc)
                        new_line = ls[0:8] + [new_desc]
                    elif cat == "transcript":
                        gene, transcript = desc.split(".")
                        # check to see if this exists in annotation
                        if annotation[gene][transcript] != None:
                            new_desc = 'transcript_id "{}.{}"; gene_id "{}"; gene_name "{}"; gene_product "{}"; gene_organism "{}";'.format(
                                gene,
                                transcript,
                                gene,
                                annotation[gene][transcript].gene,
                                annotation[gene][transcript].product,
                                annotation[gene][transcript].organism
                            )
                            if annotation[gene][transcript].db_xref:
                                db_xref = "; ".join(['ncbi_gene_identifier_{} "{}"'.format(i, j) for i,j in enumerate(annotation[gene][transcript].db_xref)])
                                new_desc = new_desc + db_xref
                        else:
                            new_desc = 'transcript_id "{}.{}"; gene_id "{}";'.format(gene, transcript, gene)
                        new_line = ls[0:8] + [new_desc]
                    elif cat in ["intron", "CDS", "exon", "start_codon", "stop_codon"]:
                        start = desc.split(";")[0]
                        name = re.search('transcript_id "(.*)"; gene_id.*', desc).groups()[0]
                        gene, transcript = name.split(".")
                        if annotation[gene][transcript] != None:
                            new_desc = 'transcript_id "{}.{}"; gene_id "{}"; gene_name "{}"; gene_product "{}"; gene_organism "{}";'.format(
                                gene,
                                transcript,
                                gene,
                                annotation[gene][transcript].gene,
                                annotation[gene][transcript].product,
                                annotation[gene][transcript].organism
                            )
                            if annotation[gene][transcript].db_xref:
                                db_xref = "; ".join(['ncbi_gene_identifier_{} "{}"'.format(i, j) for i,j in enumerate(annotation[gene][transcript].db_xref)])
                                new_desc = new_desc + db_xref
                        else:
                            new_desc = 'transcript_id "{}.{}"; gene_id "{}";'.format(gene, transcript, gene)
                        new_line = ls[0:8] + [new_desc]
                    outfile.write("{}\n".format("\t".join(new_line)))
                #pdb.set_trace()
                        

                    


if __name__ == '__main__':
	main()
