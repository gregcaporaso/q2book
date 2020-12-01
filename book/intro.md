# About this book

This book will teach you to use QIIME 2 in your microbiome research, and, if you're so inclined, to begin developing your own microbiome bioinformatics tools with QIIME 2. So grab a cup of coffee and let's get started.

## Using this book
Section 1, _Microbiome Bioinformatics with QIIME 2_, is intended to be read sequentially. This section will teach you everything you need to know to use QIIME 2 for your own microbiome research, using example data sets along the way. As you learn to use the tools, you'll be introduced to underlying theory in microbiome bioinformatics. I also sprikle in tips on how to get more out of QIIME 2 so you'll be a QIIME 2 power user by the time you're done!  

If you'd like to learn more about the underlying algorithms used in microbiome research, those are covered in Section 2 chapters. The Section 2 chapters can be ready sequentially, and derive from my earlier online book, _An Introduction to Applied Bioinformatics_, but have been re-worked to focus on microbiome-relevant examples. If you're not interested in algorithmic details, but only on learning how to use QIIME 2, you can safely skip Section 2. 

Section 3 is intended to be read sequentially. This will take you through different aspects of developing with QIIME 2. Most developers will be interested in developing plugins, so this is covered earlier in the section. We also cover developing QIIME 2 interfaces - this content will be very helpful if you'd like to integrate QIIME 2 as a component in another system you're developing. Interface development is also helpful if you have ideas about you could make QIIME 2 more accessible - dive in, and create your own graphical or other interfaces! 


## Outline (this is just a temporary holding place)

1. Microbiome Bioinformatics with QIIME 2
 a. Getting started: ways to use QIIME 2; installing QIIME 2 (refer to online documentation); download some data; do something very simple with the command line. 
 b. A very brief first tutorial: a much-abbreviated version of a web-based overview tutorial, just allowing the user to do a few simple things; view provenance
 c. A first look at core concepts of QIIME 2: semantic types, data provenance, Actions (methods, visualizers, pipelines), Results (artifacts, visualizations), importing (briefly); q2cli; artifact API (redo something from the very brief tutorial with the Artifact API).  
 d. Importing: raw sequence data; other data. 
 e. Quality control with DADA2, Deblur, others? Introduction to the FeatureTable and ASVs. 
 f. To cluster, or not to cluster? Why would or wouldn't you want to do this? Historical perspective on why we started doing this; discussion of the types of clustering supported in QIIME 2. 
 g. Taxonomy assignment: approaches; limits of resolution; training custom reference databases; non-16S data (probably the first time this matters?); taxonomy barplots; why statistics is hard on this (and we'll come back to how to do stats in the differential abundance testing chapter)
 h. Alignment and phylogenetic reconstruction: why do we care about this? what is and isn't this tree good for?
 i. Working with feature tables: filtering, ...
 j. Count normalization: rarefaction and other methods
 k. Alpha diversity: common metrics defined with worked examples; non-phylogenetic versus phylogenetic metrics; statistics; visualizations, including alpha rarefaction plots. 
 l. Beta diversity: common metrics defined with worked examples; statistics; visualizations, including PCoA plots. 
 m. Differential abundance testing
 n. Longitudinal analysis
 o. Supervised classification
 p. Network analysis? 
 q. Analysis of other -omics data types with QIIME 2: what you can do now; what you'll be able to do in the future
 r. Integrated -omics analysis

2. Understanding the algorithms (this content comes from IAB)
 a. Pairwise sequence alignment
 b. Multiple sequence alignment
 c. Sequence database searching
 d. Machine learning with Naive Bayes and Random Forest
 e. Clustering methods

3. Developing with QIIME 2
 a. Architecture of QIIME 2
 b. Building a simple plugin: qiime2.plugin
 c. Building a more complex plugin: defining types, formats, and transformers
 d. Building a simple interface: qiime2.sdk
