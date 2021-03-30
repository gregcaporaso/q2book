(alpha-diversity)=
# Alpha diversity

In this chapter we'll begin to explore metrics of microbiome diversity. We'll start with metrics of **alpha diversity**, which are measures of "within-sample" diversity. The way I typically think of these is as metrics that can be computed on a single sample.

The first sub-category of alpha diversity metric that we'll look at will be *richness*. Richness refers to how many different types of organisms are present in a sample. For example, if we're interested in species richness of plants in the Sonoran Desert and the Costa Rican rainforest, we could go to each, count the number of different species of plants that we observe, and have a basic measure of species richness in each environment.

```{admonition} Jargon: "type of organism"
As defined earlier, a "type of organism" or a "type of microbe" is an arbitrary taxonomic grouping, such as genus, species, strain, or even an amplicon sequence variant.
```

The next sub-cateogry of alpha diversity metric that we'll discuss in this chapter will be *evenness*. Evenness tells us how consistent the distribution of species abundances are in a given environment. If, for example, the most abundant plant in the Sonoran desert was roughly as common as the least abundant plant (not the case!), we would say that the evenness of plant species was high. On the other hand, if the most abundant plant was thousands of times more common than the least common plant (probably closer to the truth), then we'd say that the evenness of plant species was low.

## Measure of microbiome richness

We'll begin by looking at two measures of community richness. Both of these are commonly used in practice. 

### Observed features

Observed features is a very simple metric that can be used to quantify microbiome diversity. With this metric, we simply count the features that are observed in a given sample. A feature is considered to have been observed in a sample if it has a frequency of greater than zero. This is a **qualitative** diversity metric meaning that each feature is treated as being either observed or not observed. This metric doesn't consider how many times a feature was observed.

Let's define a feature table for this analysis:

```python
>>> sample_ids = ['A', 'B', 'C']
>>> feature_ids = ['B1','B2','B3','B4','B5','A1','E2']
>>> data = np.array([[1, 1, 5],
...                  [1, 2, 0],
...                  [3, 1, 0],
...                  [0, 2, 0],
...                  [0, 0, 0],
...                  [0, 0, 3],
...                  [0, 0, 1]])
...
>>> table2 = pd.DataFrame(data, index=feature_ids, columns=sample_ids)
>>> table2
    A  B  C
B1  1  1  5
B2  1  2  0
B3  3  1  0
B4  0  2  0
B5  0  0  0
A1  0  0  3
E2  0  0  1
```

Our sample $A$ has an observed feature frequency value of 3, sample $B$ has an observed feature frequency of 4, and sample $C$ has an observed feature frequency of 3. Note that this is different than the total counts for each column (which would be 5, 6, and 9 respectively). Based on the observed features metric, we could consider samples $A$ and $C$ to have even feature richness, and sample $B$ to have 33% higher feature richness.

We could compute this in python as follows:

```python
>>> def observed_features(table, sample_id):
...     return sum([e > 0 for e in table[sample_id]])
```

```python
>>> print(observed_features(table2, 'A'))
3
```

```python
>>> print(observed_features(table2, 'B'))
4
```

```python
>>> print(observed_features(table2, 'C'))
3
```

#### A limitation of feature counting

Imagine that we have the same table, but some additional information about the features in the table. Specifically, we've computed the following phylogenetic tree. And, for the sake of illustration, imagine that we've also assigned taxonomy to each of the features and found that our samples contain representatives from the archaea, bacteria, and eukaryotes (their labels begin with `A`, `B`, and `E`, respectively).

First, let's define a phylogenetic tree using the Newick format (which is described [here](http://evolution.genetics.washington.edu/phylip/newicktree.html), and more formally defined [here](http://evolution.genetics.washington.edu/phylip/newick_doc.html)). We'll then load that up using [scikit-bio](http://scikit-bio.org)'s [TreeNode](http://scikit-bio.org/generated/skbio.core.tree.TreeNode.html#skbio.core.tree.TreeNode) object, and visualize it with [ete3](http://etetoolkit.org).

```python
>>> import ete3
>>> ts = ete3.TreeStyle()
>>> ts.show_leaf_name = True
>>> ts.scale = 250
>>> ts.branch_vertical_margin = 15
>>> ts.show_branch_length = True
```

```python
>>> from io import StringIO
>>> newick_tree = StringIO('((B1:0.2,B2:0.3):0.3,((B3:0.5,B4:0.3):0.2,B5:0.9):0.3,'
...                        '((A1:0.2,A2:0.3):0.3,'
...                        ' (E1:0.3,E2:0.4):0.7):0.55);')
...
>>> from skbio.tree import TreeNode
...
>>> tree = TreeNode.read(newick_tree)
>>> tree = tree.root_at_midpoint()
```

```python
>>> t = ete3.Tree.from_skbio(tree, map_attributes=["value"])
>>> t.render("%%inline", tree_style=ts)
<IPython.core.display.Image object>
```

Pairing this with the table we defined above (displayed again in the cell below), given what you now know about these features, which would you consider the most diverse? Are you happy with the $\alpha$ diversity conclusion that you obtained when computing the number of observed features in each sample?

```python
>>> table2
    A  B  C
B1  1  1  5
B2  1  2  0
B3  3  1  0
B4  0  2  0
B5  0  0  0
A1  0  0  3
E2  0  0  1
```

### Phylogenetic Diversity (PD)

Phylogenetic Diversity (PD) is a metric that was developed by Dan Faith in the early 1990s (find the original paper [here](http://www.sciencedirect.com/science/article/pii/0006320792912013)). Like many of the measures that are used in microbial community ecology, it wasn't initially designed for studying microbial communities, but rather communities of "macro-organisms" (macrobes?). Some of these metrics, including PD, do translate well to microbial community analysis, while some don't translate as well. (For an illustration of the effect of sequencing error on PD, where it is handled well, versus its effect on the Chao1 metric, where it is handled less well, see Figure 1 of [Reeder and Knight (2010)](http://www.nature.com/nmeth/journal/v7/n9/full/nmeth0910-668b.html)).

PD is relatively simple to calculate. It is computed simply as the sum of the branch length in a phylogenetic tree that is "covered" or represented in a given sample. Let's look at an example to see how this works.

I'll now define a couple of functions that we'll use to compute PD.

```python
>>> def get_observed_nodes(tree, table, sample_id, verbose=False):
...     observed_features = [obs_id for obs_id in table.index
...                 if table[sample_id][obs_id] > 0]
...     observed_nodes = set()
...     # iterate over the observed features
...     for feature in observed_features:
...         t = tree.find(feature)
...         observed_nodes.add(t)
...         if verbose:
...             print(t.name, t.length, end=' ')
...         for internal_node in t.ancestors():
...             if internal_node.length is None:
...                 # we've hit the root
...                 if verbose:
...                     print('')
...             else:
...                 if verbose and internal_node not in observed_nodes:
...                     print(internal_node.length, end=' ')
...                 observed_nodes.add(internal_node)
...     return observed_nodes
...
>>> def phylogenetic_diversity(tree, table, sample_id, verbose=False):
...     observed_nodes = get_observed_nodes(tree, table, sample_id, verbose=verbose)
...     result = sum(o.length for o in observed_nodes)
...     return result
```

And then apply those to compute the PD of our three samples. For each computation, we're also printing out the branch lengths of the branches that are observed *for the first time* when looking at a given feature. When computing PD, we include the length of each branch only one time.

```python
>>> pd_A = phylogenetic_diversity(tree, table2, 'A', verbose=True)
>>> print(pd_A)
B1 0.2 0.3 0.2250000000000001 
B2 0.3 
B3 0.5 0.2 0.3 
2.0250000000000004
```

```python
>>> pd_B = phylogenetic_diversity(tree, table2, 'B', verbose=True)
>>> print(pd_B)
B1 0.2 0.3 0.2250000000000001 
B2 0.3 
B3 0.5 0.2 0.3 
B4 0.3 
2.325
```

```python
>>> pd_C = phylogenetic_diversity(tree, table2, 'C', verbose=True)
>>> print(pd_C)
B1 0.2 0.3 0.2250000000000001 
A1 0.2 0.3 0.32499999999999996 
E2 0.4 0.7 
2.65
```

How does this result compare to what we observed above with the Observed features metric? Based on your knowledge of biology, which do you think is a better representation of the relative diversities of these samples?

### Even sampling

Imagine again that we're going out to count plants in the Sonoran Desert and the Costa Rican rainforest. We're interested in getting an idea of the plant richness in each environment. In the Sonoran Desert, we survey a square kilometer area, and count 150 species of plants. In the rainforest, we survey a square meter, and count 15 species of plants. So, clearly the plant species richness in the Sonoran Desert is higher, right? What's wrong with this comparison?

The problem is that we've expended a lot more sampling effort in the desert than we did in the rainforest, so it shouldn't be surprising that we observed more species there. If we expended the same effort in the rainforest, we'd probably observe a lot more than 15 or 150 plant species, and we'd have a more sound comparison.

In sequencing-based studies of microorganism richness, the analog of sampling area is sequencing depth. If we collect 100 sequences from one sample, and 10,000 sequences from another sample, we can't directly compare the number of observed features or the phylogenetic diversity of these because we expended a lot more sampling effort on the sample with 10,000 sequences than on the sample with 100 sequences. The way this is typically handled is by randomly subsampling sequences from the sample with more sequences until the sequencing depth is equal to that in the sample with fewer sequences. If we randomly select 100 sequences at random from the sample with 10,000 sequences, and compute the alpha diversity based on that random subsample, we'll have a better idea of the relative alpha diversities of the two samples.

```python
>>> sample_ids = ['A', 'B', 'C']
>>> feature_ids = ['feature1', 'feature2', 'feature3', 'feature4', 'feature5']
>>> data = np.array([[50, 4, 0],
...                  [35, 200, 0],
...                  [100, 2, 1],
...                  [15, 400, 1],
...                  [0, 40, 1]])
...
>>> bad_table = pd.DataFrame(data, index=feature_ids, columns=sample_ids)
>>> bad_table
        A    B  C
feature1   50    4  0
feature2   35  200  0
feature3  100    2  1
feature4   15  400  1
feature5    0   40  1
```

```python
>>> print(observed_features(bad_table, 'A'))
4
```

```python
>>> print(observed_features(bad_table, 'B'))
5
```

```python
>>> print(observed_features(bad_table, 'C'))
3
```

```python
>>> print(bad_table.sum())
A    200
B    646
C      3
dtype: int64
```

### Rarefaction

#### Even sampling

#### Alpha rarefaction