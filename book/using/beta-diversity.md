(beta-diversity)=
# Beta diversity

In this chapter we'll continue our exploration of metrics of microbiome diversity. We'll next discuss metrics of **betaa diversity**, which are measures of "between-sample" diversity. I think of these as metrics that are computed based on pairs of samples. 

_content from IAB follows_

*beta diversity*) refers to **between sample diversity**, and is typically used to answer questions of the form: is sample $A$ more similar in composition to sample $B$ or sample $C$? In this section we'll explore two (of tens or hundreds) of metrics for computing pairwise dissimilarity of samples to estimate $\beta$ diversity.

### Distance metrics <link src='eac92f'/>

#### Bray-Curtis <link src='dac934'/>

The first metric that we'll look at is a quantitative non-phylogenetic $\beta$ diversity metric called Bray-Curtis. The Bray-Curtis dissimilarity between a pair of samples, $j$ and $k$, is defined as follows:

$BC_{jk} = \frac{ \sum_{i} | X_{ij} - X_{ik}|} {\sum_{i} (X_{ij} + X_{ik})}$

$i$ : feature (e.g., OTUs)

$X_{ij}$ : frequency of feature $i$ in sample $j$

$X_{ik}$ : frequency of feature $i$ in sample $k$

This could be implemented in python as follows:

```python
>>> def bray_curtis_distance(table, sample1_id, sample2_id):
...     numerator = 0
...     denominator = 0
...     sample1_counts = table[sample1_id]
...     sample2_counts = table[sample2_id]
...     for sample1_count, sample2_count in zip(sample1_counts, sample2_counts):
...         numerator += abs(sample1_count - sample2_count)
...         denominator += sample1_count + sample2_count
...     return numerator / denominator
```

```python
>>> table1
      A  B  C
OTU1  1  0  0
OTU2  3  2  0
OTU3  0  0  6
OTU4  1  4  2
OTU5  0  4  1
```

Let's now apply this to some pairs of samples:

```python
>>> print(bray_curtis_distance(table1, 'A', 'B'))
0.6
```

```python
>>> print(bray_curtis_distance(table1, 'A', 'C'))
0.857142857143
```

```python
>>> print(bray_curtis_distance(table1, 'B', 'C'))
0.684210526316
```

```python
>>> print(bray_curtis_distance(table1, 'A', 'A'))
0.0
```

```python
>>> print(bray_curtis_distance(table1, 'C', 'B'))
0.684210526316
```

Ultimately, we likely want to apply this to all pairs of samples to get a distance matrix containing all pairwise distances. Let's define a function for that, and then compute all pairwise Bray-Curtis distances between samples `A`, `B` and `C`.

```python
>>> from skbio.stats.distance import DistanceMatrix
>>> from numpy import zeros
...
>>> def table_to_distances(table, pairwise_distance_fn):
...     sample_ids = table.columns
...     num_samples = len(sample_ids)
...     data = zeros((num_samples, num_samples))
...     for i, sample1_id in enumerate(sample_ids):
...         for j, sample2_id in enumerate(sample_ids[:i]):
...             data[i,j] = data[j,i] = pairwise_distance_fn(table, sample1_id, sample2_id)
...     return DistanceMatrix(data, sample_ids)
```

```python
>>> bc_dm = table_to_distances(table1, bray_curtis_distance)
>>> print(bc_dm)
3x3 distance matrix
IDs:
'A', 'B', 'C'
Data:
[[ 0.          0.6         0.85714286]
 [ 0.6         0.          0.68421053]
 [ 0.85714286  0.68421053  0.        ]]
```

#### Unweighted UniFrac <link src='2682aa'/>

Just as phylogenetic alpha diversity metrics can be more informative than non-phylogenetic alpha diversity metrics, phylogenetic beta diversity metrics offer advantages over non-phylogenetic metrics such as Bray-Curtis. The most widely applied phylogenetic beta diversity metric as of this writing is unweighted UniFrac. UniFrac was initially presented in [Lozupone and Knight, 2005, Applied and Environmental Microbiology](http://aem.asm.org/content/71/12/8228.abstract), and has been widely applied in microbial ecology since (and the illustration of UniFrac computation presented below is derived from a similar example originally developed by Lozupone and Knight).

The unweighted UniFrac distance between a pair of samples `A` and `B` is defined as follows:

$U_{AB} = \frac{unique}{observed}$

where:

$unique$ : the unique branch length, or branch length that only leads to OTU(s) observed in sample $A$ or sample $B$

$observed$ : the total branch length observed in either sample $A$ or sample $B$

<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_d0.png" align=right/></div>

To illustrate how UniFrac distances are computed, before we get into actually computing them, let's look at a few examples. In these examples, imagine that we're determining the pairwise UniFrac distance between two samples: a red sample, and a blue sample. If a red box appears next to an OTU, that indicates that it's observed in the red sample; if a blue box appears next to the OTU, that indicates that it's observed in the blue sample; if a red and blue box appears next to the OTU, that indicates that the OTU is present in both samples; and if no box is presented next to the OTU, that indicates that it's present in neither sample.

To compute the UniFrac distance between a pair of samples, we need to know the sum of the branch length that was observed in either sample (the *observed* branch length), and the sum of the branch length that was observed only in a single sample (the *unique* branch length). In these examples, we color all of the *observed* branch length. Branch length that is unique to the red sample is red, branch length that is unique to the blue sample is blue, and branch length that is observed in both samples is purple. Unobserved branch length is black (as is the vertical branches, as those don't contribute to branch length - they are purely for visual presentation).

In the tree on the right, all of the OTUs that are observed in either sample are observed in both samples. As a result, all of the observed branch length is purple. The unique branch length in this case is zero, so **we have a UniFrac distance of 0 between the red and blue samples**.

<hr>

<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_d1.png" align=right/></div>

On the other end of the spectrum, in the second tree, all of the OTUs in the tree are observed either in the red sample, or in the blue sample. All of the observed branch length in the tree is either red or blue, meaning that if you follow a branch out to the tips, you will observe only red or blue samples. In this case the unique branch length is equal to the observed branch length, so **we have a UniFrac distance of 1 between the red and blue samples**.

<hr>

<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_d0.5.png" align=right/></div>

Finally, most of the time we're somewhere in the middle. In this tree, some of our branch length is unique, and some is not. For example, OTU 1 is only observed in our red sample, so the terminal branch leading to OTU 1 is red (i.e., unique to the red sample). OTU 2 is only observed in our blue sample, so the terminal branch leading to OTU 2 is blue (i.e., unique to the blue sample). However, the internal branch leading to the node connecting OTU 1 and OTU 2 leads to OTUs observed in both the red and blue samples (i.e., OTU 1 and OTU 2), so is purple (i.e, observed branch length, but not unique branch length). In this case, **we have an intermediate UniFrac distance between the red and blue samples, maybe somewhere around 0.5**.

<hr>
<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_with_distances.png" align=right/></div>

Let's now compute the Unweighted UniFrac distances between some samples. Imagine we have the following tree, paired with our table below (printed below, for quick reference).

```python
>>> table1
      A  B  C
OTU1  1  0  0
OTU2  3  2  0
OTU3  0  0  6
OTU4  1  4  2
OTU5  0  4  1
```

<div style="float: right; margin-left: 30px;"><img title="Image by @gregcaporaso." style="float: right; margin-left: 30px;" src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/unifrac_tree_with_distances_ab.png" align=right/></div>

First, let's compute the unweighted UniFrac distance between samples $A$ and $B$. The *unweighted* in *unweighted UniFrac* means that this is a qualitative diversity metric, meaning that we don't care about the abundances of the OTUs, only whether they are present in a given sample ($frequency > 0$) or not present ($frequency = 0$).

Start at the top right branch in the tree, and for each branch, determine if the branch is observed, and if so, if it is also unique. If it is observed then you add its length to your observed branch length. If it is observed and unique, then you also add its length to your unique branch length.

For samples $A$ and $B$, I get the following (in the tree on the right, red branches are those observed in $A$, blue branches are those observed in $B$, and purple are observed in both):

$unique_{ab} = 0.5 + 0.75 = 1.25$

$observed_{ab} = 0.5 + 0.5 + 0.5 + 1.0 + 1.25 + 0.75 + 0.75 = 5.25$

$uu_{ab} = \frac{unique_{ab}}{observed_{ab}} = \frac{1.25}{5.25} = 0.238$

As an exercise, now compute the UniFrac distances between samples $B$ and $C$, and samples $A$ and $C$, using the above table and tree. When I do this, I get the following distance matrix.

```python
>>> ids = ['A', 'B', 'C']
>>> d = [[0.00, 0.24, 0.52],
...       [0.24, 0.00, 0.35],
...       [0.52, 0.35, 0.00]]
>>> print(DistanceMatrix(d, ids))
3x3 distance matrix
IDs:
'A', 'B', 'C'
Data:
[[ 0.    0.24  0.52]
 [ 0.24  0.    0.35]
 [ 0.52  0.35  0.  ]]
```

 **TODO**: Interface change so this code can be used with ``table_to_distances``.

```python
>>> ## This is untested!! I'm not certain that it's exactly right, just a quick test.
...
... newick_tree1 = StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0),(OTU4:0.75,OTU5:0.75):1.25))root;')
>>> tree1 = TreeNode.read(newick_tree1)
...
>>> def unweighted_unifrac(tree, table, sample_id1, sample_id2, verbose=False):
...     observed_nodes1 = get_observed_nodes(tree, table, sample_id1, verbose=verbose)
...     observed_nodes2 = get_observed_nodes(tree, table, sample_id2, verbose=verbose)
...     observed_branch_length = sum(o.length for o in observed_nodes1 | observed_nodes2)
...     shared_branch_length = sum(o.length for o in observed_nodes1 & observed_nodes2)
...     unique_branch_length = observed_branch_length - shared_branch_length
...     unweighted_unifrac = unique_branch_length / observed_branch_length
...     return unweighted_unifrac
...
>>> print(unweighted_unifrac(tree1, table1, 'A', 'B'))
>>> print(unweighted_unifrac(tree1, table1, 'A', 'C'))
>>> print(unweighted_unifrac(tree1, table1, 'B', 'C'))
0.23809523809523808
0.52
0.34782608695652173
```

#### Even sampling <link src='200e13'/>

**TODO**: Add discussion on necessity of even sampling

### Interpreting distance matrices <link src='2be688'/>

In the previous section we computed distance matrices that contained the pairwise distances between a few samples. You can look at those distance matrices and get a pretty good feeling for what the patterns are. For example, what are the most similar samples? What are the most dissimilar samples?

What if instead of three samples though, we had more. Here's a screenshot from a distance matrix containing data on 105 samples (this is just the first few rows and columns):

<img src='https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/applications/images/example_big_dm.png', width=800>

Do you have a good feeling for the patterns here? What are the most similar samples? What are the most dissimilar samples?

Chances are, you can't just squint at that table and understand what's going on (but if you can, I'm hiring!). The problem is exacerbated by the fact that in modern microbial ecology studies we may have thousands or tens of thousands of samples, not "just" hundreds as in the table above. We need tools to help us take these raw distances and convert them into something that we can interpret. In this section we'll look at some techniques, one of which we've covered previously, that will help us interpret large distance matrices.

<hr>

One excellent paper that includes a comparison of several different strategies for interpreting beta diversity results is [Costello *et al.* Science (2009) Bacterial Community Variation in Human Body Habitats Across Space and Time](https://www.sciencemag.org/content/326/5960/1694.full). In this study, the authors collected microbiome samples from 7 human subjects at about 25 sites on their bodies, at four different points in time.

Figure 1 shows several different approaches for comparing the resulting UniFrac distance matrix (this image is linked from the *Science* journal website - copyright belongs to *Science*):

<img src="https://www.sciencemag.org/content/326/5960/1694/F1.large.jpg" width=800>

Let's generate a small distance matrix representing just a few of these body sites, and figure out how we'd generate and interpret each of these visualizations. The values in the distance matrix below are a subset of the unweighted UniFrac distance matrix representing two samples each from three body sites from the Costello *et al.* (2009) study.

```python
>>> sample_ids = ['A', 'B', 'C', 'D', 'E', 'F']
>>> _columns = ['body site', 'individual']
>>> _md = [['gut', 'subject 1'],
...        ['gut', 'subject 2'],
...        ['tongue', 'subject 1'],
...        ['tongue', 'subject 2'],
...        ['skin', 'subject 1'],
...        ['skin', 'subject 2']]
...
>>> human_microbiome_sample_md = pd.DataFrame(_md, index=sample_ids, columns=_columns)
>>> human_microbiome_sample_md
  body site individual
A       gut  subject 1
B       gut  subject 2
C    tongue  subject 1
D    tongue  subject 2
E      skin  subject 1
F      skin  subject 2
```

```python
>>> dm_data = np.array([[0.00, 0.35, 0.83, 0.83, 0.90, 0.90],
...                     [0.35, 0.00, 0.86, 0.85, 0.92, 0.91],
...                     [0.83, 0.86, 0.00, 0.25, 0.88, 0.87],
...                     [0.83, 0.85, 0.25, 0.00, 0.88, 0.88],
...                     [0.90, 0.92, 0.88, 0.88, 0.00, 0.50],
...                     [0.90, 0.91, 0.87, 0.88, 0.50, 0.00]])
...
>>> human_microbiome_dm = DistanceMatrix(dm_data, sample_ids)
>>> print(human_microbiome_dm)
6x6 distance matrix
IDs:
'A', 'B', 'C', 'D', 'E', 'F'
Data:
[[ 0.    0.35  0.83  0.83  0.9   0.9 ]
 [ 0.35  0.    0.86  0.85  0.92  0.91]
 [ 0.83  0.86  0.    0.25  0.88  0.87]
 [ 0.83  0.85  0.25  0.    0.88  0.88]
 [ 0.9   0.92  0.88  0.88  0.    0.5 ]
 [ 0.9   0.91  0.87  0.88  0.5   0.  ]]
```

#### Distribution plots and comparisons <link src='8fcf92'/>

First, let's look at the analysis presented in panels E and F. Instead of generating bar plots here, we'll generate box plots as these are more informative (i.e., they provide a more detailed summary of the distribution being investigated). One important thing to notice here is the central role that the sample metadata plays in the visualization. If we just had our sample ids (i.e., letters ``A`` through ``F``) we wouldn't be able to group distances into *within* and *between* sample type categories, and we therefore couldn't perform the comparisons we're interested in.

```python
>>> def within_between_category_distributions(dm, md, md_category):
...     within_category_distances = []
...     between_category_distances = []
...     for i, sample_id1 in enumerate(dm.ids):
...         sample_md1 = md[md_category][sample_id1]
...         for sample_id2 in dm.ids[:i]:
...             sample_md2 = md[md_category][sample_id2]
...             if sample_md1 == sample_md2:
...                 within_category_distances.append(dm[sample_id1, sample_id2])
...             else:
...                 between_category_distances.append(dm[sample_id1, sample_id2])
...     return within_category_distances, between_category_distances
```

```python
>>> within_category_distances, between_category_distances = within_between_category_distributions(human_microbiome_dm, human_microbiome_sample_md, "body site")
>>> print(within_category_distances)
>>> print(between_category_distances)
[0.34999999999999998, 0.25, 0.5]
[0.82999999999999996, 0.85999999999999999, 0.82999999999999996, 0.84999999999999998, 0.90000000000000002, 0.92000000000000004, 0.88, 0.88, 0.90000000000000002, 0.91000000000000003, 0.87, 0.88]
```

```python
>>> import seaborn as sns
>>> ax = sns.boxplot(data=[within_category_distances, between_category_distances])
>>> ax.set_xticklabels(['same body habitat', 'different body habitat'])
>>> ax.set_ylabel('Unweighted UniFrac Distance')
>>> _ = ax.set_ylim(0.0, 1.0)
<Figure size 432x288 with 1 Axes>
```

```python
>>> from skbio.stats.distance import anosim
>>> anosim(human_microbiome_dm, human_microbiome_sample_md, 'body site')
method name               ANOSIM
test statistic name            R
sample size                    6
number of groups               3
test statistic                 1
p-value                    0.065
number of permutations       999
Name: ANOSIM results, dtype: object
```

If we run through these same steps, but base our analysis on a different metadata category where we don't expect to see any significant clustering, you can see that we no longer get a significant result.

```python
>>> within_category_distances, between_category_distances = within_between_category_distributions(human_microbiome_dm, human_microbiome_sample_md, "individual")
>>> print(within_category_distances)
>>> print(between_category_distances)
[0.82999999999999996, 0.84999999999999998, 0.90000000000000002, 0.88, 0.91000000000000003, 0.88]
[0.34999999999999998, 0.85999999999999999, 0.82999999999999996, 0.25, 0.92000000000000004, 0.88, 0.90000000000000002, 0.87, 0.5]
```

```python
>>> ax = sns.boxplot(data=[within_category_distances, between_category_distances])
>>> ax.set_xticklabels(['same person', 'different person'])
>>> ax.set_ylabel('Unweighted UniFrac Distance')
>>> _ = ax.set_ylim(0.0, 1.0)
<Figure size 432x288 with 1 Axes>
```

```python
>>> anosim(human_microbiome_dm, human_microbiome_sample_md, 'individual')
method name                 ANOSIM
test statistic name              R
sample size                      6
number of groups                 2
test statistic           -0.333333
p-value                      0.869
number of permutations         999
Name: ANOSIM results, dtype: object
```

Why do you think the distribution of distances between people has a greater range than the distribution of distances within people in this particular example?

Here we used ANOSIM testing whether our with and between category groups differ. This test is specifically designed for distance matrices, and it accounts for the fact that the values are not independent of one another. For example, if one of our samples was very different from all of the others, all of the distances associated with that sample would be large. It's very important to choose the appropriate statistical test to use. One free resource for helping you do that is [*The Guide to Statistical Analysis in Microbial Ecology (GUSTAME)*](http://mb3is.megx.net/gustame). If you're getting started in microbial ecology, I recommend spending some time studying GUSTAME.

#### Hierarchical clustering <link src='09f456'/>

Next, let's look at a hierarchical clustering analysis, similar to that presented in panel G above. Here I'm applying the UPGMA functionality implemented in [SciPy](http://www.scipy.org/scipylib/index.html) to generate a tree which we visualize with a dendrogram. However the tips in this tree don't represent sequences or OTUs, like they did when we [covered UPGMA for phylogenetic reconstruction](alias://73d028), but instead they represent samples, and samples with a smaller branch length between them are more similar in composition than samples with a longer branch length between them. (Remember that only horizontal branch length is counted - vertical branch length is just to aid in the organization of the dendrogram.)

```python
>>> from scipy.cluster.hierarchy import average, dendrogram
>>> lm = average(human_microbiome_dm.condensed_form())
>>> d = dendrogram(lm, labels=human_microbiome_dm.ids, orientation='right',
...                link_color_func=lambda x: 'black')
<Figure size 432x288 with 1 Axes>
```

Again, we can see how the data really only becomes interpretable in the context of metadata:

```python
>>> labels = [human_microbiome_sample_md['body site'][sid] for sid in sample_ids]
>>> d = dendrogram(lm, labels=labels, orientation='right',
...                link_color_func=lambda x: 'black')
<Figure size 432x288 with 1 Axes>
```

```python
>>> labels = [human_microbiome_sample_md['individual'][sid] for sid in sample_ids]
>>> d = dendrogram(lm, labels=labels, orientation='right',
...                link_color_func=lambda x: 'black')
<Figure size 432x288 with 1 Axes>
```

### Ordination <link src='b1cdbe'/>

Finally, let's look at ordination, similar to that presented in panels A-D. The basic idea behind ordination is dimensionality reduction: we want to take high-dimensionality data (a distance matrix) and represent that in a few (usually two or three) dimensions. As humans, we're very bad at interpreting high dimensionality data directly: with ordination, we can take an $n$-dimensional data set (e.g., a distance matrix of shape $n \times n$, representing the distances between $n$ biological samples) and reduce that to a 2-dimensional scatter plot similar to that presented in panels A-D above.

Ordination is a technique that is widely applied in ecology and in bioinformatics, but the math behind some of the methods such as *Principal Coordinates Analysis* is fairly complex, and as a result I've found that these methods are a black box for a lot of people. Possibly the most simple ordination technique is one called Polar Ordination. Polar Ordination is not widely applied because it has some inconvenient features, but I find that it is useful for introducing the idea behind ordination. Here we'll work through a simple implementation of ordination to illustrate the process, which will help us to interpret ordination plots. In practice, you will use existing software, such as [scikit-bio](http://scikit-bio.org)'s [ordination module](http://scikit-bio.org/maths.stats.ordination.html).

An excellent site for learning more about ordination is [Michael W. Palmer's Ordination Methods page](http://ordination.okstate.edu/).

#### Polar ordination <link src='538e18'/>

First, let's print our distance matrix again so we have it nearby.

```python