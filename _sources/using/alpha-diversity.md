---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(alpha-diversity)=
# Alpha diversity

In this chapter we'll begin to explore metrics of microbiome diversity. We'll start with metrics of **alpha diversity**, which are measures of "within-sample" diversity. The way I typically think of these is as metrics that can be computed on a single sample.

The first sub-category of alpha diversity metric that we'll look at will be *richness*. Richness refers to how many different types of organisms are present in a sample. For example, if we're interested in species richness of plants in the Sonoran Desert and the Costa Rican rainforest, we could go to each, count the number of different species of plants that we observe, and have a basic measure of species richness in each environment.

```{admonition} Jargon: "type of organism"
As defined earlier, a "type of organism" or a "type of microbe" is an arbitrary taxonomic grouping, such as genus, species, strain, or even an amplicon sequence variant.
```

The next sub-cateogry of alpha diversity metric that we'll discuss in this chapter will be *evenness*. Evenness tells us how consistent the distribution of species abundances are in a given environment. If, for example, the most abundant plant in the Sonoran desert was roughly as common as the least abundant plant (not the case!), we would say that the evenness of plant species was high. On the other hand, if the most abundant plant was thousands of times more common than the least common plant (probably closer to the truth), then we'd say that the evenness of plant species was low.

## Measures of microbiome richness

We'll begin by looking at two measures of community richness. Both of these are commonly used in practice. 

### Observed features

Observed features is a very simple metric that can be used to quantify microbiome diversity. With this metric, we simply count the features that are observed in a given sample. A feature is considered to have been observed in a sample if it has a frequency of greater than zero. This is a **qualitative** diversity metric meaning that each feature is treated as being either observed or not observed. This metric doesn't consider how many times a feature was observed.

Let's define a simple feature table for this analysis:

```{code-cell}
:tags: [hide-input]
import pandas as pd
import numpy as np

sample_ids = ['4ac2', 'e375', '4gd8']
feature_ids = ['B1','B2','B3','B4','B5','A1','E2']
data = np.array([[1, 1, 3, 0, 0, 0, 0],
                 [1, 5, 1, 2, 0, 0, 1],
                 [3, 0, 0, 0, 0, 3, 1]])

feature_table_1 = pd.DataFrame(data, index=sample_ids, columns=feature_ids)
feature_table_1.style
```

To compute the value of observed features for each sample in our feature table, we would simply count the number of non-zero counts for each sample. All of the non-zero counts are bolded in the following view of the table. 

```{code-cell}
:tags: [hide-input]

def non_zero_blue(frequency):
    if frequency > 0:
        color = 'blue'
    else:
        color = 'black'
    return 'color: %s' % color

feature_table_1.style.applymap(non_zero_blue)
```

Counting non-zero values for each sample would result in the following values of observed features for each sample:

```{code-cell}
:tags: [hide-input]
import qiime2
import qiime2.plugins.diversity as div

feature_table_1a = qiime2.Artifact.import_data("FeatureTable[Frequency]", feature_table_1)

observed_features_1a = div.actions.alpha(feature_table_1a, metric='observed_features').alpha_diversity
observed_features_1 = observed_features_1a.view(pd.Series).to_frame(name='observed-features')
observed_features_1.style
```

Based on the observed features metric, we could consider samples `4ac2` and `4gd8` to have equal feature richness, and sample `e375` to have a higher feature richness. Also note that while `4ac2` and `4gd8` have the same feature richness, different features are present in the two samples. Richness tells us only about how many different features are present, but nothing about which features are present.

Simple counting of features, as we're doing here, is common but there are a few limitations to be aware of. 

### Even sampling

Imagine again that we're going on a sampling trip to count plants in the Sonoran Desert and the Costa Rican rainforest. We're interested in getting an idea of the plant richness in each environment. In the Sonoran Desert, we survey a square kilometer area, and count 150 species of plants. In the rainforest, we survey a square meter, and count 15 species of plants. On the basis of our survey we decide that plant species richness is higher in the desert than in the rainforest. Where did we go wrong?

We expended a lot more sampling effort in the desert than we did in the rainforest, so it shouldn't be surprising that we observed more species there. If we expended the same effort in the rainforest, we'd probably observe a lot more than 15 or 150 plant species, and we'd have a more sound comparison.

In sequencing-based studies of microbiomes, the analog of sampling area is sequencing depth. If we collect 100 sequences from one sample, and 10,000 sequences from another sample, we can't directly compare the number of observed features across these samples because we expended a lot more sampling effort on the sample with 10,000 sequences than on the sample with 100 sequences. 

The ideal way to normalize tables for computation of these metrics is a subject of ongoing research, and most likely differs depending on what you want to do with the feature table. We'll cover this topic in {ref}`feature-table-normalization`. At present, when computing alpha and {ref}`beta diversity <beta-diversity>` metics, the way this is typically handled is by randomly subsampling sequences without replacement at a fixed total frequency across all samples. This process is referred to as **rarefaction**. Continuing from the example above, if we randomly select 100 sequences to represent the sample with 10,000 sequences (i.e., we rarefy that sample to a depth of 100 sequences) we can compute its richness based on that random subsample. That richness value will serve as a more relevant comparison to the richness value for the sample that we only obtained 100 sequences from. 

```{warning}
Rarefaction is not ideal. I think of it as a necessary evil to faciliate these comparisons. Our field needs to move toward improved solutions that are easily accessible to users to get beyond the need to rarefy data. Until we have that, you **must** rarefy your data before computing alpha and beta diversity unless the metric you're using specifically does not require equal sampling effort across samples. 
```

Because rarefaction involves taking a _random_ subsample of sequences from each sample, rarefying the same feature table multiple times will yield different rarefied feature tables. This is sometimes managed by computing multiple rarefied feature tables, computing diversity metrics on each table, and then averaging the diversity value computed for each sample. We'll explore this below in {ref}`alpha-rarefaction`. The following are three different rarefied versions of our example feature table from above. Notice that the total frequency for each sample is now four - our rarefaction depth.

````{margin}
```{note}
When feature tables are rarefied in QIIME 2, features that are not observed in any samples are dropped from the table. This doesn't impact downstream analysis and results in tables with smaller file size. You can observe this in the rarefied feature tables presented in this section.
```
````

```{code-cell}
:tags: [hide-input]
import qiime2.plugins.feature_table as ft

rarefied_feature_table_1a = ft.actions.rarefy(table=feature_table_1a, sampling_depth=4).rarefied_table
rarefied_feature_table1 = rarefied_feature_table_1a.view(pd.DataFrame).astype(int)
rarefied_feature_table1.style
```

```{code-cell}
:tags: [hide-input]
rarefied_feature_table_2a = ft.actions.rarefy(table=feature_table_1a, sampling_depth=4).rarefied_table
rarefied_feature_table2 = rarefied_feature_table_2a.view(pd.DataFrame).astype(int)
rarefied_feature_table2.style
```

```{code-cell}
:tags: [hide-input]
rarefied_feature_table_3a = ft.actions.rarefy(table=feature_table_1a, sampling_depth=4).rarefied_table
rarefied_feature_table3 = rarefied_feature_table_3a.view(pd.DataFrame).astype(int)
rarefied_feature_table3.style
```

If instead of choosing to rarefy at four sequences per sample (which is lower than the total frequency of any of our samples) we rarefied at six sequences per sample (which is higher than the total frequency of sample `4ac2`), sample `4ac2` will be dropped from the resulting feature table. 

```{code-cell}
:tags: [hide-input]
rarefied_feature_table_4a = ft.actions.rarefy(table=feature_table_1a, sampling_depth=6).rarefied_table
rarefied_feature_table4 = rarefied_feature_table_4a.view(pd.DataFrame).astype(int)
rarefied_feature_table4.style
```

#### A limitation of feature counting

Imagine that we have the same table, but some additional information about the features in the table. Specifically, we've computed the following phylogenetic tree. And, for the sake of illustration, imagine that we've also assigned taxonomy to each of the features and found that our samples contain representatives from the archaea, bacteria, and eukaryotes (their labels begin with `A`, `B`, and `E`, respectively).

First, let's define a phylogenetic tree using the Newick format (which is described [here](http://evolution.genetics.washington.edu/phylip/newicktree.html), and more formally defined [here](http://evolution.genetics.washington.edu/phylip/newick_doc.html)). We'll then load that up using [scikit-bio](http://scikit-bio.org)'s [TreeNode](http://scikit-bio.org/generated/skbio.core.tree.TreeNode.html#skbio.core.tree.TreeNode) object, and visualize it with [ete3](http://etetoolkit.org).

```
import ete3
ts = ete3.TreeStyle()
ts.show_leaf_name = True
ts.scale = 250
ts.branch_vertical_margin = 15
ts.show_branch_length = True
```

```
from io import StringIO
newick_tree = StringIO('((B1:0.2,B2:0.3):0.3,((B3:0.5,B4:0.3):0.2,B5:0.9):0.3,'
                      '((A1:0.2,A2:0.3):0.3,'
                      ' (E1:0.3,E2:0.4):0.7):0.55);')

from skbio.tree import TreeNode

tree = TreeNode.read(newick_tree)
tree = tree.root_at_midpoint()
```

```
t = ete3.Tree.from_skbio(tree, map_attributes=["value"])
t.render("%%inline", tree_style=ts)

```

Pairing this with the table we defined above (displayed again in the cell below), given what you now know about these features, which would you consider the most diverse? Are you happy with the $\alpha$ diversity conclusion that you obtained when computing the number of observed features in each sample?

```
table2
```

### Phylogenetic Diversity (PD)

Phylogenetic Diversity (PD) is a metric that was developed by Dan Faith in the early 1990s (find the original paper [here](http://www.sciencedirect.com/science/article/pii/0006320792912013)). Like many of the measures that are used in microbial community ecology, it wasn't initially designed for studying microbial communities, but rather communities of "macro-organisms" (macrobes?). Some of these metrics, including PD, do translate well to microbial community analysis, while some don't translate as well. (For an illustration of the effect of sequencing error on PD, where it is handled well, versus its effect on the Chao1 metric, where it is handled less well, see Figure 1 of [Reeder and Knight (2010)](http://www.nature.com/nmeth/journal/v7/n9/full/nmeth0910-668b.html)).

PD is relatively simple to calculate. It is computed simply as the sum of the branch length in a phylogenetic tree that is "covered" or represented in a given sample. Let's look at an example to see how this works.

I'll now define a couple of functions that we'll use to compute PD.

```
def get_observed_nodes(tree, table, sample_id, verbose=False):
   observed_features = [obs_id for obs_id in table.index
               if table[sample_id][obs_id] > 0]
   observed_nodes = set()
   # iterate over the observed features
   for feature in observed_features:
       t = tree.find(feature)
       observed_nodes.add(t)
       if verbose:
           print(t.name, t.length, end=' ')
       for internal_node in t.ancestors():
           if internal_node.length is None:
               # we've hit the root
               if verbose:
                   print('')
           else:
               if verbose and internal_node not in observed_nodes:
                   print(internal_node.length, end=' ')
               observed_nodes.add(internal_node)
   return observed_nodes

def phylogenetic_diversity(tree, table, sample_id, verbose=False):
   observed_nodes = get_observed_nodes(tree, table, sample_id, verbose=verbose)
   result = sum(o.length for o in observed_nodes)
   return result
```

And then apply those to compute the PD of our three samples. For each computation, we're also printing out the branch lengths of the branches that are observed *for the first time* when looking at a given feature. When computing PD, we include the length of each branch only one time.

```
pd_A = phylogenetic_diversity(tree, table2, 'A', verbose=True)
print(pd_A)

```

```
pd_B = phylogenetic_diversity(tree, table2, 'B', verbose=True)
print(pd_B)

```

```
pd_C = phylogenetic_diversity(tree, table2, 'C', verbose=True)
print(pd_C)

```

How does this result compare to what we observed above with the Observed features metric? Based on your knowledge of biology, which do you think is a better representation of the relative diversities of these samples?



### Rarefaction

#### Even sampling

(alpha-rarefaction)=
#### Alpha rarefaction