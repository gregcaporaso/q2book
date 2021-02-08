# Getting started with using QIIME 2

This chapter will briefly introduce a few concepts that you should understand as you start using QIIME 2. This chapter is purposefully brief so you can get to the fun stuff soon. After reading this chapter you'll understand some key QIIME 2 terminology that will be used throughout this book.

## Deploying QIIME 2

The first thing you may be wondering as you prepare to learn QIIME 2 is where and how you'll deploy it. QIIME 2 can be deployed on your personal computer (e.g., your laptop or desktop computer), a cluster computer such as one owned and maintained by your university, or on cloud computing resources such as the Amazon Web Services (AWS) Elastic Compute Cloud (EC2). In the following sections I describe these options for deploying QIIME 2, and I conclude with linking you to specific instructions on how to install QIIME 2. I recommend having a working deployment of QIIME 2 as you read this book so you can run the examples yourself.

### Using QIIME 2 on your personal computer

Using QIIME 2 on your personal computer is a very convenient option, but may be too slow for some steps of your analysis workflow. Steps such as sequence quality control and taxonomic assignment may require more CPU or memory resources than are available on your personal computer. This might make some of those steps impossible to run, or might render your personal computer useless for other tasks while you're waiting for a run to complete. That could take days or even weeks, depending on how powerful your personal computer is and how big the data set is that you're working with.

A related option is that you could have a dedicated computer for running QIIME 2 analyses. This could be a server that your team owns and maintains, and that everyone on your team has access to. Then, if you need to let a job run for a long time, that machine can run for days or weeks without disrupting other work that you need to do on your personal computer. If you go with this option, it's best to think about having that computer connected to a backup power supply in case of a power outage days into an analysis. 

QIIME 2 can be installed natively on macOS or Linux personal computers, and on Windows personal computers that support the Windows Subsystem for Linux. QIIME 2 can also be used on personal computers through virtual machines. Virtual machines are a useful option if QIIME 2 can't be installed natively on your computer, but should probably not be used otherwise as the overhead of running the virtual machine will reduce the CPU and memory resources that are available to QIIME 2. 

### Using QIIME 2 on a cluster computer

If you have access to a cluster computer, this is a great option for running QIIME 2. A cluster will typically have the resources needed to run QIIME 2 on large data sets. Typically you can reach out to your institution's high performance computing or research computing office, and let them know that you need to use QIIME 2. They will have a process for having it installed on the system and made available for you to use. 

One downside of using QIIME 2 on a cluster is that often you'll only have command line (i.e., terminal) access to the cluster - not a graphical interface. This means that to view QIIME 2 results you'll need to move them off that computer to your local machine for viewing. This isn't a problem - just a minor inconvenience that you'll need to get used to.  

### Using QIIME 2 on the Amazon cloud

If you do not have access to a cluster computer, AWS is a great option for running QIIME 2. With AWS, you rent computer resources from Amazon at an hourly rate. This is very similar to running on a cluster computer, except that Amazon owns and maintains the hardware (rather than you or your institution owning and maintaining it). If you expect to run QIIME 2 analyses fairly infrequently (e.g., monthly or less) this can be a very cost effective option. You may also qualify for [grants of cloud resource time from Amazon](https://aws.amazon.com/grants/). I received a couple of generous grants from AWS when I first started exploring it for running QIIME 2.

---

```{admonition} How I run QIIME 2
:class: tip
I personally find it convenient to use both a cluster computer and my personal computer for running QIIME 2. My university has a cluster computer that all researchers at the university have access to, and our high performance computing team keeps an up-to-date version of QIIME 2 installed on that computer. I typically run the long-running steps of my workflows on the cluster and download the results to my personal computer. Then, when I'm at a more iterative stage of my analysis (for example, when generating visualizations and running statistical tests), I'll run those commands locally so I can easily view the results. Most personal computers are powerful enough for these steps.  
```

### Installing QIIME 2

Because installation instructions for QIIME 2 will change slightly between versions, these instructions are better suited for the QIIME 2 website than for a book. To get the latest installation instructions, see the [QIIME 2 installation instructions web page](https://docs.qiime2.org/2020.8/install/native/). If you have trouble installing QIIME 2, or have questions related to where or how to deploy QIIME 2, post to the [QIIME 2 Forum](https://forum.qiime2.org).

I recommend taking some time now to get QIIME 2 installed on your personal computer, or finding another way to run it so you can follow along with the examples in this book. The data sets that we'll use as we work through the examples in this book are designed to run quickly on a personal computer, so you won't need cluster or cloud access to learn from this book. If you're preparing to run a bigger analysis soon, this is probably a good time to start figuring out what computer you'll use for that work as it can take some time to work with your high performance computing office to get QIIME 2 installed on their cluster, or to get approval from your institution to bill AWS expenses to them. 

## QIIME 2 archives

One of the first things that new QIIME 2 users often notice is the `.qza` and `.qzv` files that QIIME 2 uses. All files generated by QIIME 2 are either `.qza` or `.qzv` files, and these are simply zip files that store your data alongside some QIIME 2-specific metadata. You can unzip these files with any typical unzip utility, such as WinZip, 7Zip, or unzip, and you don't need to have QIIME 2 installed to do that. 

For example, let's download a `.qza` file and take a quick look. 

```bash
$ wget https://some/url/to/sequences.qza
$ unzip sequences.qza
...
```

There are a few things I want you to notice from this example. First, we didn't use any QIIME 2 commands. We just unzipped the file as we would unzip any zip file. Next, you'll see there are two top-level directories in the output: `data` and `metadata`. The `data` directory contains a single file, `sequences.fasta`, which contains (you guessed it!) sequence data in fasta format.  If, for example, you're interested in getting some sequence data out of QIIME 2 to analyze it with another program, you can unzip your `.qza` file, and use the `sequences.fasta` file for what ever you need to do with it. Allowing users to access their data _without_ QIIME 2 was one of the earliest design goals of the system. This ensures that if QIIME 2 isn't available to you for some reason, you can still access any data that you generated with QIIME 2. You might wonder why we bothered with having QIIME 2 create these zip files in the first place, rather than just have it use the typical file formats like fasta, newick, biom, and so on. That has to do with the other information stored in the zip file - specifically the many files contained in the `metadata` directory. The files in the `metadata` directory are not intended to be viewed by a human, and you don't need any of them to work with the file (or files) in the `data` directory. QIIME 2 uses the information in the `metadata` directory to record data provenance, helping you to ensure that your bioinformatics work will be reproducible, to store a unique identifier for the data which faciliates data management, and to record the semantic type of the data (we'll come back to that shortly). This metadata empowers your bioinformatics work in ways that will be more clear as you get further along in your learning. We'll revisit this topic throughout the book. 

The `.qza` file extension is an abbreviation for QIIME Zipped Artifact, and the `.qzv` file extension is an abbreviation for QIIME Zipped Visualization. `.qza` files (which are often simply referred to as _artifacts_) are intermediary files in a QIIME 2 analysis, usually containing raw data of some sort. These files are generated by QIIME 2 and are intended to be consummed by QIIME 2. `.qzv`` files (which are often simply referred to as _visualizations_) are terminal results in a QIIME 2 analysis, such as an interactive figure or the results of a statistical test. These files are generated by QIIME 2 and are intended to be consummed by a human.

QIIME 2 also provides some of its own utilities for getting data out of `.qza` and `.qzv` files. If you're working with the QIIME 2 command line interface (which we'll use a lot in this book), the most relevant command is `qiime tools export`. If you were to run this on the `.qza` file we downloaded above, you'd see the following:

```bash
$ qiime tools export --input-path sequences.qza --output-path exported-sequences/
```

This command will unzip the archive, and take all of the files from the `data` directory and place them in `exported-sequences`. Thus if you do have QIIME 2 installed, you can get your data out of a QIIME 2 artifact without all of the QIIME 2-specific metadata using this command. 

```{admonition} Confused by the term "artifact"?
:class: note
It has been brought to our attention that the term _artifact_ can be confusing, since it is often used in science to indicate a feature that is not present naturally in a system but rather observed as a result of some technical aspect of studying that system. For example, homopolyer runs such as the `A`s in `ACTGTACTAAAAAAAAAAATGCACGTGAC` were commonly reported by some early sequencing instruments to be longer then they were in nature due to the way the sequencing reaction worked. In QIIME 2 we use the defintion of an artifact as an object that was created by some process, like an archaeological artifact. This is common in data science, and we didn't realize the potential for confusion until we were a little too far along to easily change the name.   
```  

## Plugins and actions

People generally think of QIIME 2 as a microbiome bioinformatics system, but the truth is that it's a lot more general purpose than that. QIIME 2 is built using a plugin architecture. There is a core system, which we call the _QIIME 2 Framework_ or just _The Framework_. This handles a lot of the behind-the-scenes work, like tracking data provenance and building `.qza` and `.qzv` files. There is no microbiome-specific functionality (or even bioinformatics-specific functionality) in the QIIME 2 Framework. All of the analysis functionality comes in the form of plugins to the framework. There's only a few things that you need to know about this right now. First, your deployment of QIIME 2 will have some collection of plugins installed. Plugins define actions, which are steps in an analysis workflow. For example, the `q2-diversity` QIIME 2 plugin defines actions including `alpha-phylogenetic` and `beta-phylogenetic` which can apply phylogenetic alpha and beta diversity metrics, respectively, to your data. If you don't have the `q2-diversity` plugin installed, you won't have access to those actions. To find out what QIIME 2 plugins you currently have installed, you can run the following command:

```
$ qiime --help
```

If you want to see what actions are defined by a plugin, you can call `--help` on that plugin. For example, to see what actions are available from the `q2-diversity` plugin, you can run the following command:

```
$ qiime diversity --help
```

You should see `alpha-phylogenetic` and `beta-phylogenetic` in that list, among other actions. You could go one step further if you'd like to learn about how to use those actions by calling help on the action. For example:

```
$ qiime diversity alpha-phylogenetic --help
```

The output that you get from that might look a little mysterious right now. When you finish Part 1 of this book, you'll understand how to read that help text and use it learn how to use QIIME 2 actions you've never used before. 

Another thing to know about plugins is that anyone can create and distribute them. For example, if a graduate student develops some new analysis functionality that they want to use with QIIME 2, that can create their own QIIME 2 plugin. If they want others to be able to use it, they can distribute that plugin. The [QIIME 2 Library](https://library.qiime2.org) is a website developed by the QIIME 2 team to help with dissemination of plugins. It's a great site to visit if you want to discover new analysis functionality. Developing and disseminating plugins is covered in Part 3 of this book.

QIIME 2 actions come in three varieties, as of this writing. Methods are a type of QIIME 2 action that generate one or more `.qza` files as output. Since `.qza` files are intermediary results in QIIME 2 (as discussed in the previous section), Methods typically represent some sort of processing step in your analysis, such as taxonomic annotation of sequences. The `alpha-phylogenetic` and `beta-phylogenetic` actions described above are Methods. Visualizers are a type of QIIME 2 action that generate one or more `.qzv` files as output. You'll remember that `.qzv` files are terminal results in QIIME 2, so these are steps that terminal in a workflow. An example of a QIIME 2 visualizer is the `beta-group-significance` action in the `q2-diversity` plugin, which runs PERMANOVA, ANOSIM, or PERMDISP on your data, and reports the result of the statistical test. The third type of action in QIIME 2 is a Pipeline, which can generate one or more `.qza` and/or `.qzv` as output. `Pipelines` are special in that they're a type of action that can call other actions. They are often used by developers to define simplify common workflows so they can be run by users in a single step. For example, the `core-metrics-phylogenetic` action in the `q2-diversity` plugin is a Pipeline that runs both `alpha-phylogenetic` and `beta-phylogenetic`, as well as several other actions, in a single command. In total it runs about 20 different actions, so it saves a lot of typing to be able to do achieve that with a single command.     

There's more to know about plugins and actions in QIIME 2, but this will get you started. 
   
## Semantic types 

The next topic that should be briefly covered before you start using QIIME 2 is the notion of types in QIIME 2. The term _type_ is 
overloaded with a few different concepts, so I'll start by talking about two ways that it's commonly used, and then introduce a third way that it's used less frequently but which is important to QIIME 2 (and which could help other systems, in my opinion). By disambiguating this concept now I think we'll avoid confusion later, and you'll be in a better place to understand QIIME 2 help text and other documentation. 

````{margin}
```{note}
[This video](https://www.youtube.com/watch?v=PUsvtJgpNtE) on the QIIME 2 YouTube channel discusses semantic types.
```
````

The three kinds of types that are used in QIIME 2 are file types, data types, and semantic types. File types (or formats) refer to what you probably think of when you hear that phrase: the format of a file used to store some data. For example, newick is a file type that is used for storing phylogenetic trees. Files are used most commonly for archiving data when it's not actively in use. Data types refer to how data is represented in a computer's memory (i.e., RAM) while it's actively in use. For example, if you are adding a root to an unrooted phylogenetic tree (a concept discussed in Part 2 of this book), you may use a tool like IQTree2. You would provide a path to the file containing the unrooted phylogenetic tree to IQTree2, and IQTree2 would load that tree into some _data structure_ in the computer's memory to work on it. The data structure or type, that IQTree2 uses internally to represent the phylogenetic tree will be a decision made by the developers of IQTree2. If it successfully completes the requested rooting operation, IQTree2 would write the new tree from an internal data type into a new newick-formatted file on the hard disk, and exit. As a software user, you shouldn't need to know or care about what data types are used internally by a program - you just care about what file types are used as input and output. Computer programmers care a lot about internal data types: choosing an appropriate one has huge impacts on the software. 

The third _type_ that is important in QIIME 2 is the semantic type of data. This is a representation of the _meaning_ of the data, which is not necessarily represented by either a file type or a data type. For example, two semantic types used in QIIME 2 are `Phylogeny[Rooted]` and `Phylogeny[Unrooted]`, which are used to represent rooted and unrooted trees, respectively. Both rooted and unrooted trees are commonly stored in newick files, and a computer program needs to parse (i.e., load data from a file into a in-memory data structure) to know if a tree is rooted or unrooted. For large trees, this can be a slow operation. There are some operations, such as rooting a tree, that only make sense to perform on unrooted trees. So, if you have a very large tree that you want to root, you may provide a newick file to a program that will perform that rooting. If you accidentally provide a rooted tree (say because you have tried rooting it with a few different approaches that you want to evaluate), it may take the program some time to parse the file (say 20 minutes) after which it may fail if it discovers that the tree is already rooted. That sort of delayed notification can be very frustrating as a user, since it's easily missed until a lot of time has passed. I often will start a long-running command on my university cluster computer just before the weekend. I'll typically check on the job for a few minutes, to make sure that it seems to be starting ok. I may then leave, in the hope that the job completes over the weekend and I'll have data to work with on Monday morning. It's very frustrating to come in Monday morning and find out that my job failed just a few minutes after I left on Friday for a reason that I could have quickly addressed had I known in time. 

```{note}
There's actually a worse outcome than a delayed error from a computer program when inappropriate input is provided. When a program fails and provides an error message to the user, whether or not that error message helps the user solve the problem, the program has failed loudly. Something went wrong, and it told the user about it. The program could instead fail quietly. This might happen if the program doesn't realize the input the user provided is in appropriate (e.g., an already rooted tree is provided to a program that roots an unrooted phylogenetic tree), and it runs the rooted tree through its algorithm, misinterprets something because it was provided with the wrong input, and generates an incorrect rooted tree as a result. Quiet failures can be very difficult or impossible for a user to detect, because it looks like everything has worked as expected. Failing quietly is thus _much_ worse than failing loudly - it could waste many hours of your time, and could even lead to you publishing invalid findings.
```

QIIME 2 semantic types help with this, because they provide information on what the data in a QIIME 2 `.qza` file means without having to parse anything in the `data` directory. All QIIME 2 artifacts have a semantic type associated with them (it's one of the pieces of information stored in the `metadata` directory), and QIIME 2 methods will describe what semantic types they take as input(s), and what semantic types they generate as output(s). For example, the `q2-phylogeny` plugin defines a method called `midpoint_root`. Call help on this method using the following command:

```bash
$ qiime phylogeny midpoint-root --help
```

You can see from the resulting help text that this method takes one input, an artifact of semantic type `Phylogeny[Unrooted]`. It also generates one output, an artifact of semantic type `Phylogeny[Rooted]`. This makes intuitive sense: an action that adds a root to a phylogenetic tree (as described in the help text for this method) takes an unrooted tree as input and generates a rooted tree as output. 

There is a many-to-many relationship between file types, data types, and semantic types. It's possible that a given semantic type could be represented on disk by different file types. That's well exemplified by the many different formats that are used to store demultiplexed sequence and sequence quality data. For example, this may be in one a few variants of the fastq format, or in the fasta/qual format. Additionally, data from multiple samples may be contained in one single file or split into per-sample files. Regardless of which of these file formats the data is stored in, QIIME 2 will assign the same semantic type (in this case, `SampleData[SequencesWithQuality]`. Similarly, the data type used in memory might differ depending on what operations are to be performed on the data, or based on the preference of the programmer. QIIME 2 use the semantic type `FeatureTable[Frequency]` to represent the idea of a feature table that contains counts of features (e.g., bacterial genera) on a per sample basis. Many different actions can be applied to `FeatureTable[Frequency]` artifacts in QIIME 2. When a plugin developer defines a new action that takes a `FeatureTable[Frequency]` as input, they can choose whether to load the table into a `pandas.DataFrame` or `biom.Table` object, which are two different data types.

```{note}
That last paragraph was a bit technical. Don't worry if you got lost in the details - just take away the idea that there is not a one-to-one relationship between file types, data types, and semantic types in QIIME 2. Each kind of type represents different information about the data.     

The motivation for creating QIIME 2's semantic type system was to avoid issues that can arise from providing inappropriate data to actions. The semantic type system also helps users and developers better understand the intent of QIIME 2 actions by assigning meaning to the input and output, and allows for the discovery of new potentially relevant QIIME 2 actions (more on this later). Throughout the next few chapters I'll point out semantic types of some inputs and outputs as we come across them. Again, there's more to know on this topic, but that learning can be deferred until its needed.      

## QIIME 2 View

QIIME 2 View is a web-based viewer for `.qza` and `.qzv` files. If you've never visited QIIME 2 View, take a minute to [go to the website](https://view.qiime2.org) now. This site allows for you to view QIIME 2 results on computers that don't have QIIME 2 installed on them, and there are a few examples that you can look at in the gallery on that page. I use QIIME 2 View daily, for a few different situations. First, if I'm running analyses on a cluster computer than doesn't provide a graphical interface, it's a convenient way to view those results without having to load QIIME 2 on another computer. If I have a copy of a `.qza` or `.qzv` file on my local computer (e.g., if I copied it over from the cluster) I can navigate to QIIME 2 View, and drag-and-drop the file on the QIIME 2 View page to look at it. Another scenario where this is helpful is if when I'm sharing interactive QIIME 2 results with someone who doesn't have QIIME 2 installed, for example a collaborator. I can send them `.qzv` files by email (or use the Dropbox sharing option on QIIME 2 View), and they can load and interact with the results on their computer. 

```{tip}
When you're using QIIME 2 View, your data isn't uploaded to a server. Rather the website acts as an application launcher that allows you to view local files. This means that you don't need to be concerned about data privacy issues with QIIME 2 View - the data never leaves your computer. 
```

## Getting started

Ok, that's enough discussion about QIIME 2 for now. It's time to start using it. Don't worry if you feel like you don't fully understand some of the technical details I covered in this chapter right now. My goal was to introduce these ideas to you here, and we'll revisit them in the remaining chapters in the book. As we contextualize these ideas by using QIIME 2, they'll become clearer. 
