# Convert a tertiary prediction from RGN into a PDB file
# Aleix Lafita - October 2019

library(argparse)
suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(seqinr))

###################### Argparse #############################

tertiary = "protein.tertiary"
fasta = "protein.fa"
pdb = "protein.pdb"

# create parser object
parser = ArgumentParser(
  description='Convert a tertiary prediction from RGN into a PDB file')

# specify our desired options 
parser$add_argument("-t", "--tertiary", default=tertiary,
                    help="Coordinates from RGN for protein structure [default \"%(default)s\"]")
parser$add_argument("-f", "--fasta", default=fasta,
                    help="Protein sequence in fasta format [default \"%(default)s\"]")
parser$add_argument("-p", "--pdb", default=pdb,
                    help="Name of the output pdb formatted coordinates [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

tertiary = args$tertiary
fasta = args$fasta
pdb = args$pdb


# Parse the tertiary coordinates
coords = read.csv(
  tertiary,
  sep = "",
  comment.char = "#",
  header = F,
  stringsAsFactors = F
)

seqlen = ncol(coords) / 3

coords.mat = t(coords)

pdb.df = as.data.frame(coords.mat) %>%
  mutate(
    atomr = "ATOM",
    atomid = seq(1, seqlen*3, 1),
    s1n = 7 - ceiling(log10(atomid + 1)),
    atomn = rep(c("  N   ", "  CA  ", "  C   "), seqlen),
    resn = "ALA",
    chainid = "A",
    resid = ceiling(atomid / 3),
    s2n = 4 - ceiling(log10(resid + 1)),
    x = round(V1/100, 3),
    y = round(V2/100, 3),
    z = round(V3/100, 3),
    sxn = abs(x),
    sxn = ceiling(log10(sxn)),
    sxn = ifelse(abs(x) <= 1, 1, sxn),
    sxn = 8 - ifelse(x < 0, sxn +1, sxn),
    syn = abs(y),
    syn = ceiling(log10(syn)),
    syn = ifelse(abs(y) <= 1, 1, syn),
    syn = 4 - ifelse(y < 0, syn +1, syn),
    szn = abs(z),
    szn = ceiling(log10(szn)),
    szn = ifelse(abs(z) <= 1, 1, szn),
    szn = 8 - ifelse(z < 0, szn +1, szn),
    occup = "  1.00",
    bfac = "  0.00",
    atomtype = rep(c("N", "C", "C"), seqlen)
  ) %>%
  rowwise() %>%
  mutate(
    pdbrec = paste0(
      atomr,
      paste0(rep(" ", s1n), collapse = ""),
      atomid,
      atomn,
      resn,
      " ",
      chainid,
      paste0(rep(" ", s2n), collapse = ""),
      resid,
      paste0(rep(" ", sxn), collapse = ""),
      sprintf("%.3f", x),
      paste0(rep(" ", syn), collapse = ""),
      sprintf("%.3f", y),
      paste0(rep(" ", szn), collapse = ""),
      sprintf("%.3f", z),
      occup,
      bfac,
      "           ",
      atomtype,
      "  "
    )
  )

write.table(
  pdb.df %>% select(pdbrec),
  pdb,
  row.names = F,
  col.names = F,
  quote = F
)


