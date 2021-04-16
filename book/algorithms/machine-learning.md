---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.8.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Machine learning in bioinformatics

In this chapter we'll begin talking about machine learning algorithms. Machine learning algorithms are used in bioinformatics for tasks where the user would like an algorithm to assist in the identification of patterns in a complex dataset. As is typically the case in this book, we'll work through implementing a few algorithms but these are not the implementations that you should use in practice. The code is written to be accessible for learning. [scikit-learn](http://scikit-learn.org/) is a popular and well-documented Python library for machine learning which many bioinformatics researchers and software developers use in their work. If you'd like to start trying some of these tools out, scikit-learn is a great place to start. 

```{warning}
Machine learning algorithms can easily be misused, either intentionally or unintentionally, to provide misleading results. This chapter will cover some guidelines for how to use these techniques, but it is only intended as a primer to introduce machine learning. It's not a detailed discussion of how machine learning algorithms should and shouldn't be used. If you want to start applying machine learning tools in your own research, I recommend moving from this chapter to the scikit-learn documentation, and their content on [Common pitfalls and recommended practices](https://scikit-learn.org/stable/common_pitfalls.html).
```

## The feature table

Machine learning algorithms generally are provided with a table of **samples** and user-defined **features** of those samples. These data are typically represented in a matrix, where samples are the rows and features are the columns. This matrix is referred to as a **feature table**, and it is central to machine learning and many subfields of bioinformatics. The terms used here are purposefully general. Samples are intended to be any unit of study, and features are attributes of those samples. Sometimes **labels** or **response variables** will be associated with the samples, in which case a different class of methods can be applied. 

scikit-learn provides a few example datasets that can be used for learning. Let's start by taking a look and one of them to get an idea of what input might look like in a machine learning task.

### The Iris dataset

The [Iris dataset](https://scikit-learn.org/stable/datasets/toy_dataset.html#iris-plants-dataset) is a classic example used in machine learning, originally published by RA Fisher {cite}`Fisher1936-tk`. This feature table describes four features of 150 specimens of Iris, a genus of flowering plant, representing three species. The feature table follows:

```{code-cell} ipython3
import sklearn.datasets
import pandas as pd

iris_dataset = sklearn.datasets.load_iris(as_frame=True)
iris_feature_table = iris_dataset.frame.drop('target', axis=1)
iris_feature_table.index.name = 'sample-id'
# map target integers onto species names
iris_labels = pd.Series(iris_dataset.target_names[iris_dataset.target], 
                        index=iris_dataset.target.index, name='species').to_frame()
iris_labels.index.name = 'sample-id'
```

```{code-cell} ipython3
iris_feature_table
```

The rows in this table represent our samples - in this case specimens of Iris. The columns represent features, or attributes of the samples. Each **sample vector** (i.e., row) will include a unique identifier for the sample which we usually call the _sample id_ (here these are simply integers), and values for each feature for that sample. Each **feature vector** (i.e., column) will similar contain an identifier for the feature, or the the _feature id_. These are often simplistic descriptions of the features, as they are in this example, but they don't need to be (integers would work fine as feature ids). The feature vector then contains the values measured for that feature in each sample.

This feature table on its own can serve as an input dataset for unsupervised learning tasks, which we'll cover first in this chapter. A goal of unsupervised learning might be to determine if there are clusters of samples that are most similar to one another. 

In addition to this feature table, the Iris dataset contains labels for each of the 150 samples indicating which species each sample belongs to:

```{code-cell} ipython3
iris_labels
```

The sample ids in this label vector must be the same as the sample ids in the feature table. The feature table and the sample labels together can be used as input data for supervised learning tasks, which we'll cover second in this chapter. A goal of supervised learning might be to develop a classifier that could report the species of an Iris if provided with values for its sepal length and width and its petal length and width (i.e., the features that the algorithm originally had access).

There are three different labels, or classes, in this dataset:

```{code-cell} ipython3
iris_labels['species'].unique()
```

## Unsupervised versus supervised learning methods

Many machine learning methods are classified at a high level as either unsupervised or supervised learning methods. 

In **unsupervised learning** we either don't have or don't use sample labels, and the algorithm therefore operates on a feature table alone. Typically the user is hoping to discover some structure in the data that can help them to understand which samples are most similar to each other based on their feature values. In this chapter we'll introduce ordination as an unsupervised learning task. Ordination is very widely used in biology - you may have already encountered ordination plots (such as PCoA or NMDS plots in some of your own work). 

In **supervised learning**, on the other hand, sample labels are used in addition to a feature table. As we saw above, the sample labels can be either discrete or continuous, and that distinction defines whether we're working on a classification or regression task, respectively. The goal of a supervised learning task is typically to have the computer develop a model that can accurate predict an unlabeled sample's label from its feature values (for example, what species does this Iris belong to, based on it's sepal and petal length and width).

## Machine learning methods applied to microbial sequence data

```{code-cell} ipython3
:tags: [hide-cell]

# This cell performs some configuration for this notebook. It's hidden by
# default because it's not relevant to the content of this chapter. You'll
# occasionally notice that I hide this type of information so it's not 
# distracting.

%pylab inline

import pandas as pd
import skbio
import numpy as np
import itertools
import collections
import random
```

In this chapter, we'll work with 16S rRNA data `as we did previously <load-qdr>`. Specifically, we'll load sequences from the Greengenes database and construct a feature table from them. We'll use this feature table in an unsupervised learning task and a supervised learning task. We'll also load labels for the sequences which we'll primarily use in a supervised learning task, but which we'll also use to aid in interpretation in an unsupervised learning task. 

Our goal with these tasks will be to explore phylum-level taxonomy of a few microbial phyla based on sequence data. In our unsupervised learning task, we'll determine if samples (i.e., sequences) coming from the same phyla appear to generally be more similar to each other than samples coming from different phyla. In our supervised learning task, we'll determine if we can develop a classifier to predict microbial phylum from an unlabeled sequence. 

Let's start by loading an equal number of sequences from five specific microbial phyla from Greengenes.

```{code-cell} ipython3
:tags: [hide-cell]

import qiime_default_reference as qdr
import skbio

def load_taxonomy_reference_database(phyla_of_interest, class_size=None, verbose=True):
    # Load the taxonomic data
    seq_data = {}
    phylum_to_seq_ids = {p: list() for p in phyla_of_interest}
    for e in open(qdr.get_reference_taxonomy()):
        seq_id, seq_tax = e.strip().split('\t')
        seq_tax = [e.strip() for e in seq_tax.split(';')]
        seq_phylum = ';'.join(seq_tax[:2])
            
        try:
            phylum_to_seq_ids[seq_phylum].append(seq_id)
            seq_data[seq_id] = [seq_phylum]
        except KeyError:
            # if seq_phylum is not in phylum_to_seq_ids (i.e., it
            # wasn't provided as a phylum of interest) skip this 
            # record
            pass

    for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta', 
                           constructor=skbio.DNA):
        seq_id = e.metadata['id']

        try:
            seq_data[seq_id].append(e)
        except KeyError:
            # if this seq_id wasn't previously identified as being from one of our
            # phyla of interest, skip this record
            pass
        
    if verbose:
        for phylum, seq_ids in phylum_to_seq_ids.items():
            print("%d sequences were loaded for phylum %s." % (len(seq_ids), phylum))
    
    if class_size is not None:
        sampled_seq_data = {}
        for p, seq_ids in phylum_to_seq_ids.items():
            if class_size > len(seq_ids):
                raise ValueError("Class size (%d) too large for phylum %s, which has only %d sequences." % 
                                 (class_size, p, len(seq_ids)))
            sampled_seq_ids = random.sample(seq_ids, k=class_size)
            phylum_to_seq_ids[p] = sampled_seq_ids
            sampled_seq_data.update({seq_id: seq_data[seq_id] for seq_id in sampled_seq_ids})
        seq_data = sampled_seq_data
        if verbose:
            print('\nAfter random sampling: ')
            for phylum, seq_ids in phylum_to_seq_ids.items():
                print(" %d sequences were retained for phylum %s." % (len(seq_ids), phylum))
        

    return seq_data, phylum_to_seq_ids
```

```{code-cell} ipython3
phyla = {'k__Archaea;p__Crenarchaeota', 
         'k__Archaea;p__Euryarchaeota',
         'k__Bacteria;p__Firmicutes', 
         'k__Bacteria;p__Cyanobacteria', 
         'k__Bacteria;p__Bacteroidetes', 
         'k__Bacteria;p__Actinobacteria'}
sequences_per_phylum = 100

seq_data, phylum_to_seq_ids = load_taxonomy_reference_database(phyla, sequences_per_phylum)
```

```{code-cell} ipython3
phylum_to_seq_ids['k__Archaea;p__Crenarchaeota'][:5]
```

```{code-cell} ipython3
seq_data[phylum_to_seq_ids['k__Archaea;p__Crenarchaeota'][0]][1]
```

```{code-cell} ipython3
sequence_labels = pd.Series({k:v[0] for k, v in seq_data.items()}, name='phylum').to_frame()
sequence_labels.index.name = 'id'
sequence_labels
```

```{code-cell} ipython3
sequences = pd.Series({seq_id : data[1] for seq_id, data in seq_data.items()}, 
                      name='sequence').to_frame()
sequences.index.name = 'id'
```

```{code-cell} ipython3
sequences = pd.DataFrame(seq_data).T
sequences[1][0]
```

(ml:define-nb-parameters)=

```{code-cell} ipython3
k = 7
```

```{code-cell} ipython3
kmer_frequencies = {seq_id : data[1].kmer_frequencies(k=k) for seq_id, data in seq_data.items()}
kmer_frequencies = pd.DataFrame(kmer_frequencies).fillna(0).astype(int).T
kmer_frequencies.index.name = 'id'
```

```{code-cell} ipython3
sequence_feature_table
```

```{code-cell} ipython3
# pick up here with Jaccard(?) distance computation between sequences
```

Next, we'll compute a table of the per-sequence kmer counts for all kmers in `W` for all sequences in our reference database. We'll also store the taxonomic identity of each of our reference sequences at our specified taxonomic level. We can store this information in a pandas `DataFrame`, and then view the first 25 rows of that table.

```{code-cell} ipython3
# compute all kmers for the specified alphabet
W = compute_W(alphabet, k)

# Define a function that returns the taxonomy at a specified level given
# a semi-colon separated taxonomic description.
# For example, providing 'k__Bacteria; p__Gemmatimonadetes; c__Gemm-1; o__; f__; g__; s__'
# as input will return 'k__Bacteria; p__Gemmatimonadetes' as output.
def get_taxon_at_level(taxon, level):
    taxon = [l.strip() for l in taxon.split(';')]
    return '; '.join(taxon[:level])

# Iterate over all of the reference sequences and compute their kmer frequencies.
per_sequence_kmer_counts = {}
sequence_labels = {}
for reference_sequence in reference_db:
    sequence_id = reference_sequence.metadata['id']
    
    taxon = get_taxon_at_level(reference_sequence.metadata['taxonomy'], taxonomic_level)
    sequence_labels[sequence_id] = taxon
    
    kmer_counts = dict.fromkeys(W, 0)
    kmer_counts.update(reference_sequence.kmer_frequencies(k=k))
    per_sequence_kmer_counts[sequence_id] = kmer_counts

feature_table = pd.DataFrame(data=per_sequence_kmer_counts).fillna(0).T
sequence_labels = pd.Series(sequence_labels, name='taxon')
```

```{code-cell} ipython3
# Display the first 25 samples in the feature table
feature_table[:25]
```

```{code-cell} ipython3
# Display the taxon labels for the first 25 samples
sequence_labels[:25].to_frame()
```

This table of kmer counts per taxon is our **feature table*. In this case, taxa are our samples and kmers are our features. The values in the table represent the number of times each kmer was observed in each taxon. 







## Unsupervised learning

We'll begin our exploration of machine learning approaches with unsupervised learning, and specifically discuss ordination methods. We'll work through ordination in two strokes. First, we'll explore an approach called **Polar Ordination**, where the math is simple but which isn't widely used in practice because it doesn't work well on large data sets. Working through this on a small data set will give you an idea of how ordination techniques can reduce the dimensionality of a data set and how to interpret the results of an ordination. Then, we'll apply an approach called **Principal Coordinates Analysis (PCoA)**. The math for PCoA is a bit more complicated than I want to get into in this book (I'm a biology teacher, after all), but we'll apply it to a large data set to explore how these techniques can be used in practice.

### Polar ordination

First, let's print our distance matrix again so we have it nearby.

```python
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

Polar ordination works in a few steps:

**Step 1.** Identify the largest distance in the distance matrix.

**Step 2.** Define a line, with the two samples contributing to that distance defining the endpoints.

**Step 3.** Compute the location of each other sample on that axis as follows:

$a = \frac{D^2 + D1^2 - D2^2}{2 \times D}$

where:

$D$ is distance between the endpoints

$D1$ is distance between the current sample and endpoint 1

$D2$ is distance between sample and endpoint 2.

**Step 4.** Find the next largest distance that could be used to define an *uncorrelated axis*. (This step can be labor-intensive to do by hand - usually you would compute all of the axes, along with correlation scores. I'll pick one for the demo, and we'll wrap up by looking at all of the axes.)

Here is what steps 2 and 3 look like in Python:

```python
>>> def compute_axis_values(dm, endpoint1, endpoint2):
...     d = dm[endpoint1, endpoint2]
...     result = {endpoint1: 0, endpoint2: d}
...     non_endpoints = set(dm.ids) - set([endpoint1, endpoint2])
...     for e in non_endpoints:
...         d1 = dm[endpoint1, e]
...         d2 = dm[endpoint2, e]
...         result[e] = (d**2 + d1**2 - d2**2) / (2 * d)
...     return d, [result[e] for e in dm.ids]
```

```python
>>> d, a1_values = compute_axis_values(human_microbiome_dm, 'B', 'E')
>>> for sid, a1_value in zip(human_microbiome_dm.ids, a1_values):
...     print(sid, a1_value)
A 0.0863586956522
B 0
C 0.441086956522
D 0.431793478261
E 0.92
F 0.774184782609
```

```python
>>> d, a2_values = compute_axis_values(human_microbiome_dm, 'D', 'E')
>>> for sid, a2_value in zip(human_microbiome_dm.ids, a2_values):
...     print(sid, a2_value)
A 0.371193181818
B 0.369602272727
C 0.0355113636364
D 0
E 0.88
F 0.737954545455
```

```python
>>> from pylab import scatter
>>> ord_plot = scatter(a1_values, a2_values, s=40)
<Figure size 432x288 with 1 Axes>
```

And again, let's look at how including metadata helps us to interpret our results.

First, we'll color the points by the body habitat that they're derived from:

```python
>>> colors = {'tongue': 'red', 'gut':'yellow', 'skin':'blue'}
>>> c = [colors[human_microbiome_sample_md['body site'][e]] for e in human_microbiome_dm.ids]
>>> ord_plot = scatter(a1_values, a2_values, s=40, c=c)
<Figure size 432x288 with 1 Axes>
```

And next we'll color the samples by the person that they're derived from. Notice that this plot and the one above are identical except for coloring. Think about how the colors (and therefore the sample metadata) help you to interpret these plots.

```python
>>> person_colors = {'subject 1': 'red', 'subject 2':'yellow'}
>>> person_c = [person_colors[human_microbiome_sample_md['individual'][e]] for e in human_microbiome_dm.ids]
>>> ord_plot = scatter(a1_values, a2_values, s=40, c=person_c)
<Figure size 432x288 with 1 Axes>
```

#### Determining the most important axes in polar ordination <link src='fb483b'/>

Generally, you would compute the polar ordination axes for all possible axes. You could then order the axes by which represent the largest differences in sample composition, and the lowest correlation with previous axes. This might look like the following:

```python
>>> from scipy.stats import spearmanr
...
>>> data = []
>>> for i, sample_id1 in enumerate(human_microbiome_dm.ids):
...     for sample_id2 in human_microbiome_dm.ids[:i]:
...         d, axis_values = compute_axis_values(human_microbiome_dm, sample_id1, sample_id2)
...         r, p = spearmanr(a1_values, axis_values)
...         data.append((d, abs(r), sample_id1, sample_id2, axis_values))
...
>>> data.sort()
>>> data.reverse()
>>> for i, e in enumerate(data):
...     print("axis %d:" % i, end=' ')
...     print("\t%1.3f\t%1.3f\t%s\t%s" % e[:4])
axis 0: 	0.920	1.000	E	B
axis 1: 	0.910	0.943	F	B
axis 2: 	0.900	0.928	E	A
axis 3: 	0.900	0.886	F	A
axis 4: 	0.880	0.543	E	D
axis 5: 	0.880	0.429	F	D
axis 6: 	0.880	0.429	E	C
axis 7: 	0.870	0.371	F	C
axis 8: 	0.860	0.543	C	B
axis 9: 	0.850	0.486	D	B
axis 10: 	0.830	0.429	C	A
axis 11: 	0.830	0.406	D	A
axis 12: 	0.500	0.232	F	E
axis 13: 	0.350	0.143	B	A
axis 14: 	0.250	0.493	D	C
```

So why do we care about axes being uncorrelated? And why do we care about explaining a lot of the variation? Let's look at a few of these plots and see how they compare to the plots above, where we compared axes 1 and 4.

```python
>>> ord_plot = scatter(data[0][4], data[1][4], s=40, c=c)
<Figure size 432x288 with 1 Axes>
```

```python
>>> ord_plot = scatter(data[0][4], data[13][4], s=40, c=c)
<Figure size 432x288 with 1 Axes>
```

```python
>>> ord_plot = scatter(data[0][4], data[14][4], s=40, c=c)
<Figure size 432x288 with 1 Axes>
```

#### Interpreting ordination plots <link src='40e0a6'/>

There are a few points that are important to keep in mind when interpreting ordination plots. Review each one of these in the context of polar ordination to figure out the reason for each.

**Directionality of the axes is not important (e.g., up/down/left/right)**

One thing that you may have notices as you computed the polar ordination above is that the method is *not symmetric*: in other words, the axis values for axis $EB$ are different than for axis $BE$. In practice though, we derive the same conclusions regardless of how we compute that axis: in this example, that samples cluster by body site.

```python
>>> d, a1_values = compute_axis_values(human_microbiome_dm, 'E', 'B')
>>> d, a2_values = compute_axis_values(human_microbiome_dm, 'E', 'D')
>>> d, alt_a1_values = compute_axis_values(human_microbiome_dm, 'B', 'E')
```

```python
>>> ord_plot = scatter(a1_values, a2_values, s=40, c=c)
<Figure size 432x288 with 1 Axes>
```

```python
>>> ord_plot = scatter(alt_a1_values, a2_values, s=40, c=c)
<Figure size 432x288 with 1 Axes>
```

Some other important features:

* Numerical scale of the axis is generally not useful
* The order of axes is generally important (first axis explains the most variation, second axis explains the second most variation, ...)
* Most techniques result in uncorrelated axes.
* Additional axes can be generated (third, fourth, ...)

### Principle Coordinates Analysis (PCoA)

```
import qiime2
import qiime2.plugins.feature_table as ft
import qiime2.plugins.diversity as div

# Iterate over all of the reference sequences and compute their kmer frequencies.
per_sequence_kmer_counts = []
for reference_sequence in reference_db:
    #taxon = get_taxon_at_level(reference_sequence.metadata['taxonomy'], taxonomic_level)
    kmer_counts = dict.fromkeys(W, 0)
    kmer_counts.update(reference_sequence.kmer_frequencies(k=k))
    per_sequence_kmer_counts.append(pd.Series(kmer_counts, name=reference_sequence.metadata['id']))

# Build a table of the kmer frequencies as a pandas.DataFrame object, and then 
# display the first 25 rows of that table.
per_sequence_kmer_counts = pd.DataFrame(data=per_sequence_kmer_counts).fillna(0).T

feature_table_1a = qiime2.Artifact.import_data("FeatureTable[Frequency]", per_sequence_kmer_counts.T)
```

```
jaccard_1a = div.actions.beta(feature_table_1a, metric='jaccard').distance_matrix
```

```
taxa_of_interest = ['k__Archaea', 'p__Cyanobacteria', 'p__Firmicutes', 'p__Bacteroidetes', 'p__Proteobacteria']
metadata = {}
for reference_sequence in reference_db:
    id_ = reference_sequence.metadata['id']
    taxon = get_taxon_at_level(reference_sequence.metadata['taxonomy'], taxonomic_level)
    label_as_other = True
    for taxon_of_interest in taxa_of_interest:
        # this approach is horrendous
        if taxon_of_interest in taxon:
            label_as_other = False
    if label_as_other:
        metadata[id_] = 'Other'
    else:
        metadata[id_] = taxon
metadata = pd.Series(metadata, name='taxon').to_frame()
metadata.index.name = 'id'
metadata = qiime2.Metadata(metadata)
```

```
pcoa_1a = div.actions.pcoa(jaccard_1a).pcoa
```

```
import qiime2.plugins.emperor as emperor

emperor.actions.plot(pcoa_1a, metadata).visualization
```

```
import skbio.stats.ordination
ordination = pcoa_1a.view(skbio.stats.ordination.OrdinationResults)

_ = ordination.plot(metadata.to_dataframe(), column='taxon', cmap='Set1')
```

## Supervised classification

We'll continue our exploration of machine learning approaches with **supervised classification**, and specifically with an algorithm called **Naive Bayes**.  We'll implement Naive Bayes to gain an understanding of how it works, and I think you'll discover that this idea of machines learning isn't quite as mysterious or science fiction-y as it sounds. The math involved in Naive Bayes is relatively straight-forward, which is why I chose this algorithm to present here. There are many machine algorithms with more complex math, but Naive Bayes is  widely used and powerful, so it's a good place to get started. 

We'll explore supervised classification in the context of a now familiar topic: taxonomic classification of 16S rRNA sequences. We previously explored this problem in {doc}`database-searching`, so it's worth spending a few minutes skimming that chapter if it's not fresh in your mind.

Briefly, the problem that we are going to address here is as follows. We have a query sequence ($q_i$) which is not taxonomically annotated (meaning we don't know the taxonomy of the organism whose genome it is found in), and a reference database ($R$) of taxonomically annotated sequences ($r_1, r_2, r_3, r_n$). We want to infer a taxonomic annotation for $q_i$. We'll again work with [Greengenes](http://greengenes.secondgenome.com/), a 16S rRNA sequence database, which we'll access using [QIIME default reference project](https://github.com/biocore/qiime-default-reference). (This should all sound very familiar - if not, I again suggest that you review {doc}`database-searching`.)

Before we get to this though, lets talk about what supervised classification algorithms are and how the classifiers they build are evaluated. 

### Defining a classification task

In a classification task, there are two or more pre-defined classes, and the goal is to assign observations to those classes. As humans, we run perform these kinds of tasks everyday. For example, if you're browsing a bookstore you might classify titles as ones you want to read versus everything else (the ones you're not interested in reading). You might group the apps that you have on your phone into folders by classifying them by category (e.g., "school", "entertainment", or "social media"). 

When we're working with large data sets, supervised classification algorithms can help us with classification tasks that will make us more efficient or help us understand our data. A classic example of this outside of bioinformatics is an email spam filter. For every email that is received, the spam filter must define it as spam or not spam so the message can directed either to the user's spam folder or the user's inbox. The stakes can be high: a filter that is too permissive will cause the user's inbox to get filled with junk mail, while a filter that is overly restrictive could cause relevant messages to be directed to the spam folder. In either case, the email user could miss important messages.

In the case of taxonomic assignment, our classes will be taxonomic groups at a user-defined taxonomic level. For example, a phylum classifier for 16S rRNA sequences would take an unannotated sequence as input and as output present the phylum that the sequence most likely originated from.

### Training data, test data, and cross-validation

Supervised classification algorithms need to be provided with data that is used to develop a model to use in classification (in other words, to train the classifier). This data is a collection of observations with defined classes, and is referred to as the **training data**. These labeled examples are the "supervision" aspect of supervised learning. In the email spam filter example, this would be email messages that are annotated as either spam or not spam. In the 16S taxonomy assignment example, this would be 16S sequences that are taxonomically annotated. It is typically important that the training data be balanced - in other words, that there are roughly the same number of examples of each class.

In addition to the training data, an independent collection of observations with defined classes is needed as **test data**. These observations are not used to train the classifier, but rather to evaluate how the classifier performs on previously unseen data. The goal of testing the classifier on these test data is to predict what performance will be on **real world** data. Real world data refers to data for which the class is currently unknown. In the spam filter example, real world data would be new emails that you are receiving. In the 16S rRNA taxonomy assignment example, real world data could be sequences that you obtain from the environment using a DNA sequencing instrument. The test data shouldn't be used for optimization of classifiers: in other words, you shouldn't develop a classifier on training data, test it on test data, go back and make changes to the classifier, and then re-test on test data. This would risk **over-fitting** the classifier to a particular test data set and performance on that test data may no longer be predictive of how the classifier will perform when it is used on real world data. 

Because training and test data sets can be very costly to develop (for example, they may require many hours of annotation by humans) we often use an approach call **k-fold cross validation** during classifier development and optimization {numref}`cross-validation-1`. In k-fold cross-validation, the training data is split into `k` different data sets, where `k` is usually five or ten. In each of the data sets, $1/k$ of the entries are used as test data and all of the other entries are used as training data. In `k` iterations, the classifier is developed on the training data and tested on the test data. The average performance of the classifier is then computed across the `k` iterations. k-fold cross validation therefore allows for developing and optimizing a classifier without using dedicated test data.

```{figure} ./images/ml-cross-validation.png
---
name: cross-validation-1
---
An illustration of k-fold cross validation where a single data set is split into k independent training and test data sets. Each circle represents a labeled entry for use in training or testing, and colors indicate the class of each entry. In the case of a spam filter, for example, red circles might represent spam messages while green circles represent messages that are not spam.
Image source: [Gufosowa](https://commons.wikimedia.org/wiki/File:K-fold_cross_validation_EN.svg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), via Wikimedia Commons.
```

### Evaluating a binary classifier 

As mentioned above, in a classification task there are two or more pre-defined classes. A binary classifier would be a specific type of classifier for which there are exactly two classes - for example, spam and not spam. We'll start talking about how classifiers are evaluated by discussing binary classifiers because they're the easiest to understand. 

Imagine we're building a classifier that attempts to predict whether an individual is healthy or has some specific disease (let's call it *Disease X*). Perhaps the data that the classifier uses is based on a variety of medical data that has undergone a feature extraction process to generate features that can be used by a supervised classification algorithm. When a classifier is developed, you can think of it like a function that will take a collection of features for a sample and return a value of "healthy" or "diseased". 

The goal of our classifier is to serve as a diagnostic tool that identifies whether a patient has Disease X based on features of their medical data. A positive test result therefore indicates that the patient has Disease X while a negative test result indicates that they are healthy. When we apply our classifier to test data (i.e., where we know the correct class), there are a few possible outcomes. 

* The classifier predicts a positive test result, and the sample is known to come from a patient with Disease X. This is a **true positive (TP)**. 
* The classifier predicts a positive test result, and the sample is known to come from a healthy patient. This is a **false positive (FP)**. FPs are also referred to as type 1 errors.
** The classifier predicts a negative test result, and the sample is known to come from a patient with Disease X. This is a **false negative (FN)**. FNs are also referred to as type 2 errors.
* The classifier predicts a negative test result, and the sample is known to come from a healthy patient. This is a **true negative (TN)**. 

A classifier would typically be evaluated by running it on many samples and tallying the count of TP, FP, FN, and TN results. These tallies are typically presented in a structure known as a **confusion matrix**. For the confusion matrix, there many different values that can be computed which inform us of some aspect of classifier performance. 

The simplest way to think about evaluating the performance of our classifier from a confusion matrix is to compute its **accuracy** as:

```{math}
:label: accuracy
accuracy = \frac{TP + TN}{TP + FP + FN + TN}
```

In words, accuracy can be defined as the fraction of the total test cases that the classifier classified correctly. Accuracy gives us an idea of the classifier performance, but it hides some potentially relevant information from us. A low accuracy classifier could, for example, almost never achieve false positives but frequently achieve false negatives. Such a classifier could still be a clinically useful tool. Because false positives are very infrequent but false negatives are common, that means when the classifier indicates a positive test result that person nearly always has the disease. If the classifier indicates a negative result, that could be an indicator that additional testing is needed. Of course we would rather our classifier achieve fewer false negatives, but if this is a very cheap test and the additional tests are more expensive, it could be a useful first screening approach. 

Two other metrics are more widely used for evaluating classifiers, and these are typically computed as a pair. These metrics are **precision** and **recall** and they are more informative than accuracy because they indicate whether a classifier might suffer more from false positives or false negatives. 

Precision is the fraction of the positives reported by the classifier that are actually positives, or:

```{math}
:label: precision
precision = \frac{TP}{TP + FP}
```

Recall is the fraction of the actual positives that are reported to be positive by the classifier, or:

```{math}
:label: recall
recall = \frac{TP}{TP + FN}
```

Precision thus tells us how frequently our classifier yields false positives, while recall tells us how frequently our classifier yields false negatives. We of course would always like both of these values to be high, but depending on the application of our classifier, we may prefer high precision over high recall, or we may prefer high recall over high precision. 

### Naive Bayes classifiers

Naive Bayes classifiers work by building a model of what different classes look like based on labeled training data. In the case of taxonomic assignment of 16S sequences, the classes are different microbial taxonomy names at a given taxonomic level. For example, Proteobacteria and Cyanobacteria would be two classes if we were building a classifier for bacterial phyla. The data that is provided would be the 16S sequences associated with different representatives of the classes, but more specifically Naive Bayes needs these sequences broken into finer-grained features for it to work well. The development of features from raw training data is referred to as **feature extraction**. This can be part of the classifier training software, or it can be independent. Features of sequences could be nearly anything, such as sequence length, presence or absence of certain sequence patterns (or motifs), GC content, and so on. The most commonly used features for sequence classification tasks such as this is {ref}`overlapping kmers <kmer>`, which we have previously seen when looking at heuristic algorithms for database searching. In this case, feature extraction for a given sequence would involve the identification of all of the kmers contained in that sequence.

In this chapter, instead of using sequence alignment to identify the most likely taxonomic origin of a sequence, we'll train Naive Bayes classifiers to do this by building {ref}`kmer <kmer>`-based models of the 16S sequences of taxa in our reference database. We'll then run our query sequences through those models to identify the most likely taxonomic origin of each query sequence. Since we know the taxonomic origin of our query sequences in this case, we can evaluate the accuracy of our classifiers by seeing how often they return the known taxonomy assignment. If our training and testing approaches are well-designed, the performance on our tests will inform us of how accurate we can expect our classifier to be on data where the actual taxonomic origin is unknown. 



### Training a Native Bayes classifier 

The first thing our Naive Bayes classifier will need is the set of all possible words of length ``k``. This will be dependent on the value of ``k`` and the characters in our alphabet (i.e., the characters that we should expect to find in the reference database). This set is referred to as ``W``, and can be computed as follows.

```{code-cell} ipython3
alphabet = skbio.DNA.nondegenerate_chars
k = 2

def compute_W(alphabet, k):
    return set(map(''.join, itertools.product(alphabet, repeat=k)))

W = compute_W(alphabet, k)

print('Alphabet contains the characters: %s' % ', '.join(alphabet))
print('For an alphabet size of %d, W contains %d length-%d kmers.' % (len(alphabet), len(W), k))
```

```{admonition} Exercise
Given the DNA alphabet (A, C, G, and T), how many different kmers of length 3 are there (i.e., 3-mers)? How many different 7-mers are there? How many 7-mers are there if there are twenty characters in our alphabet (as would be the case if we were working with protein sequences instead of DNA sequences)?
```

The next thing we'll need to train our classifier is a way to extract all kmers from a given sequence. scikit-bio provides this functionality in the ``skbio.DNA`` sequence object (as well as in the other sequence object types). It also provides functionality for computing the kmer frequencies in a given sequence. This information can be obtained for one of our reference sequences as follows:

```{code-cell} ipython3
kmers = reference_db[0].iter_kmers(k=k)
for kmer in kmers:
    print(kmer, end=' ')
```

That's a lot of kmers, and of course many of them are present multiple times. Tallies of the frequencies of each kmer can be computed as follows.

```{code-cell} ipython3
print(reference_db[0].kmer_frequencies(k=k))
```

This information can be convenient to store in a pandas ``Series`` object:

```{code-cell} ipython3
pd.Series(reference_db[0].kmer_frequencies(k=k), name=reference_db[0].metadata['id'])
```

To train our taxonomic classifier, we next need to define a few things. First, at what level of taxonomic specificity do we want to classify our sequences? We should expect to achieve higher accuracy at less specific taxonomic levels such as phylum or class, but these are likely to be less informative biologically than more specific levels such as genus or species. Let's start classifying at the phylum level to keep our task simple, since we're working with a small subset of the reference database here. In Greengenes, phylum is the second level of the taxonomy.

Next, how long should our kmers be? We don't have a good idea of this to start with. The longer our kmers, the more likely they are to be specific to certain taxa, which is good because that will help with classification. However, if they get too long it becomes less likely that we'll observe those kmers in sequences that aren't represented in our database because the longer the sequence is the more likely we are to see variation across other organisms that are assigned to the same taxonomy. Based on some of my own work in this area, I'll start us out with 7-mers (i.e., kmers of length 7).

Finally, we'll need to know the value of `W`, defined above as the set of all possible kmers given our alphabet and the value of `k`.



With this information, we'll next compute our kmer probability table. The content of this table will be the probability of observing every kmer in W given a taxon. This is computed based on a few values:

$N$ : The total number of sequences in the training set (i.e., our reference database).

$W$: The set of all possible kmers, given $k$ and an alphabet.

$w_i$: An individual kmer in $W$.

$n(w_i)$ : The number of total sequences containing $w_i$.

$P_i$ : The probability of observing $w_i$. Initially it might seem as though this would be computed as $n(w_i) / N$, but this neglects the possibility that a kmer observed in a query sequence might not be represented in our reference database (i.e., $n(w_i) = 0$), which would create problems later, when we're assigning probabilities to each class for query sequences. As a result, 0.5 is added to the numerator and 1 is added to the denominator. When we alter counts in this way, we refer to the values that we're adding as **pseudocounts**. 

$P(w_i | taxon)$ : The probability of observing a kmer given a taxon. Again, it would seem that this would be computed as the proportion of sequences in the taxon containing the kmer, but this would neglect that we'll likely observe kmers in our query sequences that are not represented in our reference database. A pseudocount is therefore added again to the numerator and denominator. This time the pseudocount in the numerator is scaled by how frequent the kmer is in the reference database as a whole: specifically, it is $P_i$.

Our "kmer probability table" is $P(w_i | taxon)$ computed for all kmers in W and all taxa represented in our reference database. We'll compute that and again look at the first 25 rows.

```{code-cell} ipython3
def compute_kmer_probability_table(feature_table, sequence_labels):
    N = feature_table.shape[0] # number of training sequences

    # number of sequences containing kmer wi
    n_wi = feature_table.astype(bool).sum(axis=0)
    n_wi.name = 'n(w_i)'

    # probabilities of observing each kmer
    Pi = (n_wi + 0.5) / (N + 1)
    Pi.name = 'P_i'
    
    # number of times each taxon appears in training set
    taxon_counts = collections.Counter(sequence_labels)

    
    taxon_table = feature_table.astype(bool).groupby(by=sequence_labels, axis=0).sum()
    
    # probabilities of observing each kmer in each taxon
    p_wi_t = []
    for taxon, count in taxon_counts.items():
        p_wi_t.append(pd.Series((taxon_table.loc[taxon] + Pi) / (count + 1), name=taxon))

    return pd.DataFrame(p_wi_t).T
```

```{code-cell} ipython3
kmer_probability_table = compute_kmer_probability_table(feature_table, sequence_labels)
kmer_probability_table[:25]
```

```{raw-cell}
This kmer probability table represents our kmer-based models of the phyla in our reference database. We can use this table to compute probabilities of taxonomically unannotated query sequences belonging to each of the phyla represented in this table.

### Applying a Naive Bayes classifier 

With our kmer probability table we are now ready to classify unknown sequences. We'll begin by defining some query sequences. We'll pull these at random from our reference sequences, which means that some of the query sequences will be represented in our reference database and some won't be. We'll also trim out 200 bases of our reference sequences since (as of this writing) we typically don't obtain full-length 16S sequences from a DNA sequencing instrument. We're thus trying to emulate a real-world classification of environmental 16S rRNA sequences, where some might be perfect matches to sequences we've observed before while others might represent previously unobserved sequences (or even previously unknown organisms).
```

```{code-cell} ipython3
:tags: [hide-cell]

def load_taxonomy_query_sequences(start_position=100, length=200):
    queries = []
    for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta', constructor=skbio.DNA):
        if e.has_degenerates():
            # For the purpose of this lesson, we're going to ignore sequences that contain
            # degenerate characters (i.e., characters other than A, C, G, or T)
            continue
        e = e[start_position:start_position + length]
        queries.append(e)

    return queries
```

We'll load a collection of query sequences as we did in {doc}`database-searching`.

```{code-cell} ipython3
import random

queries = load_taxonomy_query_sequences()
queries = random.sample(queries, k=50)
```

Again, we can index into these results to look at individual sequences. Note that because we're trying to emulate working with unannotated sequences here, the query sequences don't have taxonomic annotations in their metadata.

```{code-cell} ipython3
queries[0]
```

For a given query sequence, its taxonomy will be classified as follows. First, the set of all kmers will be extracted from the sequence. This is referred to as $V$. Then, for all taxa in the kmer probability table, the probability of observing the query sequence will be computed given that taxon: $P(query | taxon)$. This is computed as the product of all its kmer probabilities for the given taxon. (It should be clear based on this formula why it was necessary to add pseudocounts when computing our kmer probability table - if not, kmer probabilities of zero would result in a zero probability of the sequence being derived from that taxon at this step.)

After computing $P(query | taxon)$ for all taxa, the taxonomy assignment returned is simply the one achieving the maximum probability. Here we'll classify a sequence and look at the resulting taxonomy assignment.

```{code-cell} ipython3
# This function classifies a sequence that has already been split into a list
# of kmers.
def classify_V(V, kmer_probability_table):
    P_S_t = [] # probability of the sequence given the taxon
    for taxon in kmer_probability_table:
        kmer_probabilities = kmer_probability_table[taxon]
        probability = 1.0
        for v_i in V:
            probability *= kmer_probabilities[v_i]
        P_S_t.append((probability, taxon))
    return max(P_S_t)[1], V

# This function is a little more convenient to use. It classifies a sequence 
# directly, first by computing V, and then by calling classify_V.
def classify_sequence(query_sequence, kmer_probability_table, k):
    V = list(map(str, query_sequence.iter_kmers(k=k)))
    return classify_V(V, kmer_probability_table)
```

```{code-cell} ipython3
query = queries[0]
taxon_assignment, V = classify_sequence(query, kmer_probability_table, k)
print(taxon_assignment)
```

Since we know the actual taxonomy assignment for this sequence, we can look that up in our reference database. Was the assignment correct?

```{code-cell} ipython3
get_taxon_at_level(reference_taxonomy[query.metadata['id']], taxonomic_level)
```

```{admonition} Exercise
Try classifying a few other query sequences and determining if the returned class was correct. You can do this by changing which entry in `queries` you're assigning to the value `query`. Keep track of how many times the classifier returned the correct assignment.
```

### Evaluating our confidence in the results of the Naive Bayes classifier

Because the query and reference sequences that were working with were randomly selected from the full reference database, each time you run this notebook you should observe different results. Chances are however that if you run the above steps multiple times you'll get the wrong taxonomy assignment at least some of the time. Up to this point, we've left out an important piece of information: how confident should we be in our assignment, or in other words, how dependent is our taxonomy assignment on our specific query? If there were slight differences in our query (e.g., because we observed a very closely related organism, such as one of the same species but a different strain, or because we sequenced a different region of the 16S sequence) would we obtain the same taxonomy assignment? If so, we should have higher confidence in our assignment. If not, we should have lower confidence in our assignment. This is additionally important because our classifier will _always_ return one of the classes, even if our query sequence is very different than any of the sequences in our reference database.

We can quantify confidence using an approach called bootstrapping. With a bootstrap approach, we'll get our taxonomy assignment as we did above, but then for some user-specified number of times, we'll create random subsets of V sampled with replacement. We'll then assign taxonomy to each random subset of V, and count the number of times the resulting taxonomy assignment is the same as the one we received when assigning taxonomy to V. The count of times that they are the same divided by the number of iterations we've chosen to run will be our confidence value. If the assignments are often the same we'll have a high confidence value. If the assignments are often different, we'll have a low confidence value.

Let's now assign taxonomy and compute a confidence for that assignment.

```{code-cell} ipython3
def classify_sequence_with_confidence(sequence, kmer_probability_table, k,
                                      confidence_iterations=100):
    # classify the query sequence, as we did above
    taxon, V = classify_sequence(sequence, kmer_probability_table, k)

    count_same_taxon = 0
    # Define the size of each subsample as 10% of the actual number of
    # kmers in the query sequence.
    subsample_size = int(len(V) * 0.1)
    # Perform n iterations (where n is provided by the user as 
    # confidence_iterations) where a random subset of the query sequence's
    # kmers are used for the classification task.
    # Keep track of the number of times the observed result is the same as
    # that for the query sequence. 
    for i in range(confidence_iterations):
        subsample_V = np.random.choice(V, subsample_size, replace=True)
        subsample_taxon, _ = classify_V(subsample_V, kmer_probability_table)
        if taxon == subsample_taxon:
            count_same_taxon += 1
    confidence = count_same_taxon / confidence_iterations

    return (taxon, confidence)
```

```{code-cell} ipython3
taxon_assignment, confidence = classify_sequence_with_confidence(queries[0], kmer_probability_table, k)
print(taxon_assignment)
print(confidence)
```

How did the computed confidence compare to the accuracy taxonomy assignment?

At first glance, we don't necessarily have an idea of what good versus bad confidence scores are, but we can use our reference database to explore that. Knowing that can allows us to develop a confidence threshold that we can use in our work. For example, we can define a confidence threshold above which we would accept a taxonomy assignment and below which we might reject it. To explore this, let's compute taxonomy assignments and confidence for all of our query sequences and then see what the distributions of confidence scores look like for correct assignments and incorrect assignments.

```{code-cell} ipython3
correct_assignment_confidences = []
incorrect_assignment_confidences = []
summary = []

for query in queries:
    predicted_taxonomy, confidence = classify_sequence_with_confidence(query, kmer_probability_table, k)
    actual_taxonomy = get_taxon_at_level(reference_taxonomy[query.metadata['id']], taxonomic_level)
    if actual_taxonomy == predicted_taxonomy:
        correct_assignment_confidences.append(confidence)
    else:
        incorrect_assignment_confidences.append(confidence)

    summary.append([predicted_taxonomy, actual_taxonomy, confidence])
summary = pd.DataFrame(summary, columns=['Predicted taxonomy', 'Actual taxonomy', 'Confidence'])
```

```{code-cell} ipython3
import seaborn as sns

ax = sns.boxplot(data=[correct_assignment_confidences, incorrect_assignment_confidences])
ax = sns.swarmplot(data=[correct_assignment_confidences, incorrect_assignment_confidences], color="black")
_ = ax.set_xticklabels(['Correct assignments', 'Incorrect assignments'])
_ = ax.set_ylabel('Confidence')

ax
```

What does this plot tell you about how well setting a confidence threshold is likely to work? If you never wanted to reject a correct assignment, how often would you accept an incorrect assignment? If you never wanted to accept an incorrect assignment, how often would you reject a correct assignment?

```{admonition} Exercise
Jump back up to where we [defined `k` and `taxonomic_level`](ml:define-nb-parameters) and modify those values. How does the accuracy of the classifier change if you increase or decrease `k` while keeping the value of `taxonomic_level` fixed? How does the accuracy change if you increase or decrease the `taxonomic_level` while keeping `k` fixed? 
```

## Variations on the input to machine learning algorithms

As in the Iris dataset, the labels in our microbial data are discrete (i.e., categorical or qualitative) as opposed to continuous (i.e., quantitative). If our labels in a supervised learning project were continous instead of discrete - for example the abundance of an organism in an environment - we could still supervised learning, but we would work with different algorithms. Specifically, we'd used supervised regression algorithms, rather than supervised classification algorithms.  

Similarly, while the features we worked with in our unsupervised and supervised learning examples were continuous values, feature values could also be discrete (e.g., the sex of a subject, or the species of a specimen in an environment). The applicable algorithms might change, but machine learning techniques in general would still be available. 

scikit-learn provides other example datasets, including [the diabetes dataset](https://scikit-learn.org/stable/modules/generated/sklearn.datasets.load_diabetes.html#sklearn.datasets.load_diabetes), [the housing market dataset](https://scikit-learn.org/stable/datasets/toy_dataset.html#boston-house-prices-dataset) and [the hand-writing dataset](https://scikit-learn.org/stable/datasets/toy_dataset.html#optical-recognition-of-handwritten-digits-dataset). These are good illustrations of other types of data that can be used in machine learning tasks. The message to take away is that if you can wrangle your data into a feature table, potentially with corresponding sample labels, you will likely be able to apply machine learning techniques to that data. That said, and as I mentioned at the beginning of this chapter, this introduction barely scratches the surface of this complex branch of statistics and computer science. Especially with the accessible of these methods through software like scikit-learn, it's easy to get to the point where you know enough to get yourself into trouble by using machine learning methods inappropriately. If you'd like to apply these tools in your research, you must continue your learning. I recommend continuing with [scikit-learn's documentation](https://scikit-learn.org/).

## List of works cited

```{bibliography} ../references.bib
:filter: docname in docnames
```

```{code-cell} ipython3

```
