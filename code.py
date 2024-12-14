!pip install pyranges pandas seaborn

# download human gtf
!wget --content-disposition -c https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.chr.gtf.gz

# download cat gtf
!wget --content-disposition -c  https://ftp.ensembl.org/pub/release-112/gtf/felis_catus/Felis_catus.Felis_catus_9.0.112.chr.gtf.gz

# download chicken gtf
!wget --content-disposition -c  https://ftp.ensembl.org/pub/release-112/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.chr.gtf.gz

# download Ehuxleyi gtf
!wget --content-disposition -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-59/gtf/emiliania_huxleyi/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.59.gtf.gz

# download yeast gtf
!wget https://ftp.ensembl.org/pub/release-112/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz

# uncompress
!gunzip *.gz

# explore human gtf
!head Homo_sapiens.GRCh38.112.chr.gtf

# explore e huxleyi gtf
!head Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.59.gtf

# imports
import pyranges as pr
import pandas as pd
import seaborn as sns
import numpy as np

# read all files as pyranges object at once

human_gr = pr.read_gtf("Homo_sapiens.GRCh38.112.chr.gtf")
cat_gr = pr.read_gtf("Felis_catus.Felis_catus_9.0.112.chr.gtf")
chicken_gr = pr.read_gtf("Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.chr.gtf")

ehuxleyi_gr = pr.read_gtf("Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.59.gtf")
yeast_gr = pr.read_gtf("Saccharomyces_cerevisiae.R64-1-1.112.gtf")

def total_genes(gr):
  """Return total genes in a gtf"""
  genes = gr[gr.Feature == "gene"]
  return len(genes)

# DO NOT DELETE

all_genes = pd.DataFrame({"species" :["human",
                                      "cat",
                                      "chicken",
                                      "yeast",
                                      "ehuxleyi"],
                          "total_genes": [total_genes(human_gr),
                                          total_genes(cat_gr),
                                          total_genes(chicken_gr),
                                          total_genes(yeast_gr),
                                          total_genes(ehuxleyi_gr)]})
all_genes# DO NOT DELETE

def genes_by_genebiotype(gr):
  """Return the number of genes by gene_biotype"""
  genes = gr[gr.Feature == "gene"]
  return genes.gene_biotype.value_counts() # DO NOT DELETE

genes_by_genebiotype(human_gr) # DO NOT DELETE

def total_protein_coding_genes(gr):
  """Get total protein coding genes in  gtf"""
  genes = gr[gr.Feature == "gene"]
  genes = genes[genes.gene_biotype == "protein_coding"]
  return len(genes) # DO NOT DELETE

# DO NOT DELETE
pc_genes = pd.DataFrame({"species" :["human", "cat", "chicken", "yeast", "ehuxleyi"],
                      "total_genes": [total_protein_coding_genes(human_gr), total_protein_coding_genes(cat_gr),
                                      total_protein_coding_genes(chicken_gr),  total_protein_coding_genes(yeast_gr), total_protein_coding_genes(ehuxleyi_gr)]})
pc_genes # DO NOT DELETE

def genes_by_genebiotype_and_strand(gr):
  """Returns number of genes by strand keeping only protein_coding, lncRNA and processed_pseudogene.
  Returns a dataframe with cols [gene_biotype, Strand, count]
  """
  genes = gr[gr.Feature == "gene"]
  df = genes.df.loc[:,['gene_biotype', 'Strand']].value_counts()
  df = df.reset_index()
  df = df.loc[df.gene_biotype.isin(['protein_coding', 'lncRNA', 'processed_pseudogene'])]
  return df

genes_by_genebiotype_and_strand(human_gr) # DO NOT DELETE

genes_by_genebiotype_and_strand(chicken_gr) # DO NOT DELETE

genes_by_genebiotype_and_strand(ehuxleyi_gr) # DO NOT DELETE

def number_of_isoforms_per_gene(gr):
  """Return the number of unique isoforms per gene
  Returns a dataframe with cols [gene_id, number_of_isoforms]
  """
  genes = gr[gr.Feature == "transcript"]
  genes = genes[genes.gene_biotype == "protein_coding"]
  df = pd.DataFrame(genes.df.gene_id.value_counts()).reset_index()
  df.columns = ["gene_id", "number_of_isoforms"]
  return df

human_nisoforms_per_gene = number_of_isoforms_per_gene(human_gr) # DO NOT DELETE
human_nisoforms_per_gene # DO NOT DELETE

# DO NOT DELETE

grs = {"human": human_gr, "cat": cat_gr, "chicken": chicken_gr, "yeast": yeast_gr, "ehuxleyi": ehuxleyi_gr}
nexons_per_isoform_list = []
for species, gr in grs.items():
  df = number_of_isoforms_per_gene(gr)
  df["species"] = species
  nexons_per_isoform_list.append(df)

nisoforms_per_gene_df = pd.concat(nexons_per_isoform_list)
sns.violinplot(data=nisoforms_per_gene_df, x="number_of_isoforms", y="species", hue="species")

# DO NOT DELETE

nisoforms_per_gene_subset1 = nisoforms_per_gene_df[nisoforms_per_gene_df.species.isin(["yeast", "ehuxleyi"])]
sns.violinplot(data=nisoforms_per_gene_subset1, x="number_of_isoforms", y="species", hue="species")

# DO NOT DELETE

nisoforms_per_gene_subset2 = nisoforms_per_gene_df[~nisoforms_per_gene_df.species.isin(["yeast", "ehuxleyi"])]
sns.violinplot(data=nisoforms_per_gene_subset2, x="number_of_isoforms", y="species", hue="species")

def get_isoform_length(gr, gene_biotype):
  """Return a dataframe with length of all isoforms of a certain gene_biotype
  Your output should be a series of lengths (i.e. nothing else is necessary)
  """
  gr = gr[gr.Feature == "transcript"]
  gr = gr[gr.gene_biotype == gene_biotype].df
  gr["Length"] = gr.End - gr.Start
  return gr["Length"]

# INSERT CODE HERE TO PLOT THE DISTRIBUTION OF length of human isoforms of type:
# protein coding, lncRNA and processed pseudogene on a violinplot
# The xaxis of violin plot should be log10(length) and xaxis should be the name
# of the transcript

protein_coding = get_isoform_length(human_gr, gene_biotype="protein_coding")
lncrna = get_isoform_length(human_gr, gene_biotype="lncRNA")
pseudogene = get_isoform_length(human_gr, gene_biotype="processed_pseudogene")
protein_coding = pd.DataFrame({"length": protein_coding, "type": ["protein_coding"] * len(protein_coding)})
lncrna = pd.DataFrame({"length": lncrna, "type": ["lncRNA"] * len(lncrna)})
pseudogene = pd.DataFrame({"length": pseudogene, "type": ["processed_pseudogene"] * len(pseudogene)})

df = pd.concat([protein_coding, lncrna, pseudogene])
df["log10_length"] = np.log10(df["length"])
sns.violinplot(data=df, x="log10_length", y="type", hue="type")


def number_of_exons_per_isoform(gr):
  """Return number of exons per isoform as a dataframe.
  Returns a dataframe with cols [transcript_id, number_of_exons]
  """
  gr = gr[gr.gene_biotype == "protein_coding"]
  gr_df = gr.df
  gr_df = gr_df.loc[:,["transcript_id", "exon_id"]].drop_duplicates()
  gr_df = gr_df.groupby("transcript_id").count().reset_index()
  gr_df.columns = ["transcript_id", "number_of_exons"]
  return gr_df

# DO NOT DELETE
grs = {"human": human_gr, "cat": cat_gr, "chicken": chicken_gr, "yeast": yeast_gr, "ehuxleyi": ehuxleyi_gr}
nexons_per_isoform_list = []
for species, gr in grs.items():
  df = number_of_exons_per_isoform(gr)
  df["species"] = species
  nexons_per_isoform_list.append(df)

nexons_per_isoform_df = pd.concat(nexons_per_isoform_list)
sns.violinplot(data=nexons_per_isoform_df, x="number_of_exons", y="species", hue="species")
del nexons_per_isoform_df
del nexons_per_isoform_list

def longest_isoform_per_gene(gr):
  """Get the longest isoform per protein coding gene.
  Returns a dataframe with one isoform per gene id along with isoform length
  """

  gr = gr[gr.gene_biotype == "protein_coding"]
  gr = gr[gr.Feature == "transcript"]
  gr_df = gr.df
  gr_df['tx_length'] = gr_df.End - gr_df.Start
  idx = gr_df.groupby("gene_id")["tx_length"].idxmax()
  gr_df = gr_df.iloc[idx]
  return gr_df

# DO NOT DELETE
grs = {"human": human_gr, "cat": cat_gr, "chicken": chicken_gr, "yeast": yeast_gr, "ehuxleyi": ehuxleyi_gr}
longest_isoform_list = []
for species, gr in grs.items():
  df = longest_isoform_per_gene(gr)
  df["species"] = species
  longest_isoform_list.append(df)

longest_isoform_list_df = pd.concat(longest_isoform_list)
longest_isoform_list_df["log10_tx_length"] = np.log10(longest_isoform_list_df["tx_length"])
sns.violinplot(data=longest_isoform_list_df, x="tx_length", y="species", hue="species")


sns.violinplot(data=longest_isoform_list_df, x="log10_tx_length", y="species", hue="species") # DO NOT DELETE


def get_longest_isoform_cds_lengths_per_gene(gr, transcript_id_list):
  """Get the longest isoform CDS length per protein coding gene
  Given a list of transcript_ids, returns a dataframe with columns gene_id and cds_length (total length of coding domain sequence of the gene)
  """

  gr_cds = gr[gr.Feature == "CDS"].df
  gr_cds = gr_cds.loc[gr_cds.transcript_id.isin(transcript_id_list)]
  gr_cds['Length'] = gr_cds.End - gr_cds.Start
  # sum length per gene_id
  gr_cds_gene = gr_cds.groupby("gene_id").Length.sum().reset_index()
  gr_cds_gene.columns = ["gene_id", "cds_length"]
  return gr_cds_gene

# DO NOT DELETE
cds_len_list = []
for species, gr in grs.items():
  tx_list = longest_isoform_list_df.loc[longest_isoform_list_df.species == species].transcript_id.tolist()
  df = get_longest_isoform_cds_lengths_per_gene(gr, tx_list)
  df["species"] = species
  cds_len_list.append(df)

cds_len_df = pd.concat(cds_len_list)
cds_len_df["cds_log10_length"] = np.log10(cds_len_df["cds_length"])
sns.violinplot(data=cds_len_df, x="cds_length", y="species", hue="species")

sns.violinplot(data=cds_len_df, x="cds_log10_length", y="species", hue="species")
