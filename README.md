**General description**

PhyloTax is an R function designed to assign taxonomy to unclassified ASVs/OTUs based on their phylogenetic relationships with closely related sequences of known taxonomy.

Traditional taxonomic assignment methods such as the Naive Bayesian Classifier (NBC) and Bayesian Lowest Common Ancestor (BLCA) often leave a proportion of ASVs/OTUs unassigned, particularly when reference databases are incomplete. PhyloTax complements these approaches by using phylogenetic information to infer taxonomy for these unclassified sequences.

**How it works**

The function identifies ASVs that are unassigned at a user-defined taxonomic level and searches for their closest classified relatives in a phylogenetic tree.

For each unassigned ASV:

The ASV is located within the phylogenetic tree.
The parent node containing the ASV is identified.
A subset of the tree is extracted using:
R
tree.sub <- tree_subset(tree, node = parent, levels_back = LEVELBACK)
Show more lines
Starting with LEVELBACK = 0, the function progressively expands the subset tree until it contains ASVs with assigned taxonomy.
The taxonomy shared among all classified relatives in the subset tree is determined.
Taxonomic assignment is then made at the lowest taxonomic rank that is consistently shared by all relatives, provided that the node supporting the relationship exceeds the specified bootstrap threshold.
Example

If all classified relatives belong to the species Nitrosomonas europaea:

Kingdom	Bacteria
Phylum	Pseudomonadota
Class	Betaproteobacteria
Order	Nitrosomonadales
Family	Nitrosomonadaceae
Genus	Nitrosomonas
Species	Nitrosomonas europaea-like

If the closest classified relatives belong to different species within the same genus (e.g. Nitrosomonas europaea and Nitrosomonas communis), the ASV is assigned as:

Kingdom	Bacteria
Phylum	Pseudomonadota
Class	Betaproteobacteria
Order	Nitrosomonadales
Family	Nitrosomonadaceae
Genus	Nitrosomonas_like
Species Unassigned Nitrosomonas_like

PhyloTax assumes that all classified relatives used for assignment have consistent taxonomic information. However, inconsistencies can occur when taxonomies originate from multiple sources or from NBC/BLCA classifications.

To address this issue, PhyloTax can:

Verify taxonomic consistency among neighbouring classified ASVs.
Use a reference taxonomy (reftax) instead of the ASV taxonomy table.
Query the NCBI Taxonomy database using the classification() function from the taxize package when inconsistencies are detected.

This helps ensure that phylogeny-based assignments are based on coherent and up-to-date taxonomy.

Function syntax

PhyloTax(
tax = ASV_taxonomy,
bootstrap_treshold = 0.75,
check_consistency = TRUE,
unassigned_level = "Phylum",
reftax = Ref.Tax,
tree = tree,
useNCBI = TRUE)

**Arguments**

_tax_: Taxonomy table containing ASV classifications. Row names must correspond to ASV names, and columns should contain taxonomic ranks (Kingdom to Species).

_reftax_ (data.frame): Reference taxonomy table used to validate and harmonise taxonomic assignments. Row names must correspond to ASV names and columns should contain: Column 1: Kingdom; Column 2: Phylum; Column 3: Class; Column 4: Order; Column 5: Family; Column 6: Genus; Column 7: Species

_tree_ (phylo object): A phylogenetic tree object containing all ASVs. Tip labels must correspond to row names in the taxonomy table.

_bootstrap_treshold_ (numeric): Minimum bootstrap support required for phylogeny-based assignment. (Default: 0.75)

_check_consistency_ (Logical): whether taxonomic consistency among classified relatives should be checked before assignment. (Default: TRUE)

_useNCBI_ (Logical): whether the NCBI Taxonomy database should be queried when inconsistencies are detected in the reference taxonomy. (Default: TRUE)

_unassigned_level_ (character): Taxonomic rank used to identify unassigned ASVs. Possible values include: "Kingdom"; "Phylum"; "Class"; "Order"; "Family"; "Genus"; "Species" (Default: "Phylum")

**Requirements**:
ape/ treeio/ tidytree/ taxize
