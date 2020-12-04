# About this book

There is rapidly growing interest in microbiome research, including research into human microbiomes (the trillions of microorganisms that live in and on our bodies and which are integral to human health), environmental microbiomes (such as those found in soil which are directly relevant to climate change and agricultural research), animal microbiomes, food microbiomes, and others. We are only beginning to glimpse the range of applications of microbiome science. 

The technologies that are enabling work in these fields are the same that are driving the data revolution in biology. Primarily this work is driven by high-throughput DNA sequencing, which is applied for profiling microbial community composition (e.g., marker gene profiling such as 16S or ITS sequencing), functional potential (such as shotgun metagenomic sequencing), and functional activity (such as metatranscriptome sequencing). Other “omics” technologies are now playing an increasing role in microbiome research, such as mass-spectrometry-based metabolomics, which provides profiles of small molecule metabolites in an environment, and metaproteomics which provides more detailed descriptions of functional activities of microbes (and their hosts, if applicable). As a result, bioinformatics software tools are essential to microbiome research. For many microbiome researchers, bioinformatics is an intimidating and challenging aspect of their projects. 

My focus for the past decade has been on developing tools to empower researchers to analyze their microbiome data. This work is centered around the QIIME microbiome bioinformatics platform, and currently QIIME 2. The original QIIME, now referred to as [QIIME 1](http://qiime.org), was published in 2010 {cite}`qiime1` and has now been cited over 22,000 times in the primary literature (Google Scholar, November 2020). QIIME 2, which was published in July of 2019 {cite}`qiime2`, has succeeded QIIME 1 and has over 1500 citations already (including citations of the pre-print). QIIME 2 is better than QIIME 1 in all ways, and QIIME 1 is no longer actively supported. If you have previously used QIIME 1, you should invest time in learning and switching to QIIME 2. If you're new to QIIME, start with QIIME 2. 

QIIME 2 has large and growing user and developer communities, and these communities make QIIME 2 possible. The epicenter of the community is the [QIIME 2 Forum](https://forum.qiime2.org). The forum is primarily known as a place where users can get technical support with QIIME 2 for no charge. Developers of QIIME 2 moderate the forum, and typically repond to technical support questions within a couple of business days. The forum is also a great place to discuss general topics in microbiome bioinformatics, or microbiome research methods generally. There are many active discussions on these topics on the forum. Keeping up with the discussions on the forum is a great way to learn about current topics in microbiome research methods. There's also a free job board on the forum - you can use the forum to find jobs, or post your own job ads there to find employees who are well-versed in QIIME 2 and other bioinformatics tools. If you're not already a member of the QIIME 2 Forum, you should consider joining. It's a great way for you to get help, and as you develop your QIIME 2 skills helping others on the forum is a great way to reenforce your learning and to get involved in the community.

The QIIME 2 project is popular, essential to microbiome research, and (as of this writing) stabley funded through federal and other grants. Time spent reading this book and learning QIIME 2 will help you to advance your microbiome research projects. This book can be used for self-learning, as a companion text in [QIIME 2 workshops](https://workshops.qiime2.org), or even as a text in university courses that cover QIIME 2. 

So grab a cup of coffee and let's get started! ☕

## Using this book
Section 1, _Microbiome Bioinformatics with QIIME 2_, is intended to be read sequentially. This section will teach you everything you need to know to get started using QIIME 2 for your own microbiome research, using example data sets along the way. As you learn to use the tools, you'll be introduced to underlying theory in microbiome bioinformatics. I also sprikle in tips and tricks on how to get more out of QIIME 2 so you'll be a QIIME 2 power user by the time you're done!  

If you'd like to learn more about the underlying algorithms used in microbiome research, some of those are covered in Section 2 chapters. The Section 2 chapters can be read sequentially but don't have to be. These chapters derive from my earlier online book, _An Introduction to Applied Bioinformatics_, but have been re-worked to focus on microbiome-relevant problems. If you're not interested in algorithmic details, but only on learning how to use QIIME 2, you can safely skip Section 2. Many of the algorithms discussed in these chapters are fundamentals of bioinformatics, so if you're interested in becoming a bioinformatician or bioinformatics software developer, you should read these.  

Section 3 is intended to be read sequentially. This will take you through different aspects of developing with QIIME 2. Most developers will be interested in creating plugins, so building, documenting, testing, and distributing QIIME 2 plugins is covered first. Building QIIME 2 plugins is a great way to get your bioinformatics tools in the hands of a lot of microbiome researchers. It comes with other benefits as well: QIIME 2's unique retrospective data provenance tracking system will be used whenever your plugins are run, ensuring reproducibility of the work and providing information that will help you provide technical support to users; your users will be able to access your functionality through any of QIIME 2's interfaces, which provide access to the same tools through interfaces geared toward users with different levels of computational sophisitcation; and QIIME 2 plugins are often published as stand-alone papers, so building QIIME 2 plugins can help you get the publications you may need to support your career. Section 3 also covers developing QIIME 2 interfaces. This content will help you integrate QIIME 2 as a component in another system you're developing, or allow you to implement your own ideas to make QIIME 2 more accessible. If you are considering wrapping QIIME 2 in another system, this content will give you tools to simplify that work through the QIIME 2 software development kit (SDK). QIIME 2 and the development team enthusiastically support your plugin and interface development efforts! Don't forget to get in touch on the forum with your QIIME 2 development questions.  


## Outline (this is just a temporary holding place)

1. Microbiome Bioinformatics with QIIME 2
	* Getting started: ways to use QIIME 2; installing QIIME 2 (refer to online documentation); download some data; do something very simple with the command line. 
	* A very brief first tutorial: a much-abbreviated version of a web-based overview tutorial, just allowing the user to do a few simple things; view provenance
	* A first look at core concepts of QIIME 2: semantic types, data provenance, Actions (methods, visualizers, pipelines), Results (artifacts, visualizations), importing (briefly); q2cli; artifact API (redo something from the very brief tutorial with the Artifact API).  
	* Importing: raw sequence data; other data. 
	* Quality control with DADA2, Deblur, others? Introduction to the FeatureTable and ASVs. 
	* To cluster, or not to cluster? Why would or wouldn't you want to do this? Historical perspective on why we started doing this; discussion of the types of clustering supported in QIIME 2. 
	* Taxonomy assignment: approaches; limits of resolution; training custom reference databases; non-16S data (probably the first time this matters?); taxonomy barplots; why statistics is hard on this (and we'll come back to how to do stats in the differential abundance testing chapter)
	* Alignment and phylogenetic reconstruction: why do we care about this? what is and isn't this tree good for?
	* Working with feature tables: filtering, ...
	* Count normalization: rarefaction and other methods
	* Alpha diversity: common metrics defined with worked examples; non-phylogenetic versus phylogenetic metrics; statistics; visualizations, including alpha rarefaction plots. 
	* Beta diversity: common metrics defined with worked examples; statistics; visualizations, including PCoA plots. 
	* Differential abundance testing
	* Longitudinal analysis
	* Supervised classification
	* Network analysis? 
	* Analysis of other -omics data types with QIIME 2: what you can do now; what you'll be able to do in the future
	* Integrated -omics analysis

2. Understanding the algorithms (this content comes from IAB)
	 * Pairwise sequence alignment
	 * Multiple sequence alignment
	 * Sequence database searching
	 * Machine learning with Naive Bayes and Random Forest
	 * Clustering methods

3. Developing with QIIME 2
	 * Architecture of QIIME 2
	 * Building a simple plugin: qiime2.plugin, usage API
	 * Building a more complex plugin: defining types, formats, and transformers
	 * Building a simple interface: qiime2.sdk
