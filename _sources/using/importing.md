(importing)=
# Importing data into QIIME 2

A QIIME 2 analysis almost always starts with importing data for use in QIIME 2. This step creates [QIIME 2 archives](getting-started:archives) from data in other formats, such as fastq or biom. To import data into QIIME 2, you need to define the [file type and semantic type](getting-started:semantic-types) of the data. 

I'll be honest with you: in addition to being the first step in a user's QIIME 2 analysis, importing is often the most challenging step. The reason it's challenging is that there are tens or even hundreds of different file types that users would like to use with QIIME 2, and many file formats in bioinformatics are poorly defined. For example, the [ete3 phylogenetic analysis and visualization toolkit](http://etetoolkit.org/) recognizes (as of this writing) [11 different variants of the newick file format](http://etetoolkit.org/docs/latest/reference/reference_tree.html#ete3.TreeNode). A newick file doesn't include information in it on which of these variants it is, so it's up to the person working with the file to know that. A user importing data into QIIME 2 needs to have a very good understanding of what format their data is in, and then learn how to provide that information to QIIME 2. 

This chapter will provide an overview of importing data into QIIME 2, focused on a few of the most commonly imported formats. After reading this content, you can find examples of importing additional data types in [](importing-examples).

## Why is importing necessary?

More directly relevant for importing data into QIIME 2 is the fastq file format. The fastq file format stores sequence and associated sequence quality information, and uses a clever approach for storing quality information that enables the quality to be represented in the same number of characters as the sequence itself. For example, here is a single sequence and quality record from a fastq file:

```
@M00176:65:000000000-A41FR:1:1101:9905:3163 1:N:0:0
AACCAGCACCTCAAGTGGTCAGGATGATTATTGGGCCTAAAGCATCCGTAGCCGGATCTGTAAGTTTTCGGTTAAATCTGTACGCTCAACGTACAGGCTGCCGGGAATACTGCAGATCTAGGGAGTGGGAGAGGTAGACGGTACTCGGTAG
+
AHAABBABFBFFGGGDGBGGGECFGHHHHHHHHGHHGGHHHHHFHHHGFHGGHGGGGGHHHHHFHHHHHGGGGGHHHHHGHHHHFGEEGHGHHHGGHGHGGHGGGGGHHHHHHHHHHHHFHHGGGCFFGHGGGGFFDGGFG<GEHHGGG/C
```

````{margin}
```{note}
If you'd like to learn more about the `fastq` format, see the [scikit-bio documentation](http://scikit-bio.org/docs/latest/generated/skbio.io.format.fastq.html) and the [Wikipedia entry](https://en.wikipedia.org/wiki/FASTQ_format).
```
````

The line beginning with the `@` symbol indicates the beginning of a new sequence record. It is followed by an identifier for this sequence that, in this example, was generated during an Illumina MiSeq sequencing run. The next line contains the sequence. The line beginning with the `+` symbol indicates the end of the sequence, and the last line indicates the quality of each base call in the sequence. Each of the characters on this line represents an encoded Phred quality score. For example, in this fastq file `A` might represent a quality score of 32, and `H` might represent a quality score of 39. You can refer to a simple translation table [such as this one](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm), to decode the quality scores. That seems simple enough - so what's the problem? Well, in another fastq file, `A` might represent a quality score of 1, and `H` might represent a quality score of 8. Again, you could look those values up in a different translation table. But the problem is that the fastq file itself doesn't contain explicit information about what encoding scheme was used, so when trying to interpret the information in the file, without additional context you won't know if `A` represents a high quality base call or a low quality base call. Ouch! There are some approaches that can be applied to infer how scores are encoded, but they are not completely reliable and it can be computationally expensive to figure out. Again, the burden is on the person working with the fastq file to know [which encoding scheme](https://en.wikipedia.org/wiki/FASTQ_format#Encoding) is used. 

One of the core design goals of QIIME 2 was that it should keep track of the meaning of data in the files it's using, like how quality scores are encoded in fastq files. This removes that burden from the user, and ensures that someone who encounters the data at a later time (for example, you or your boss in five years) will know how to interpret it. Continuing with the example of fastq files, because the quality score encoding scheme isn't stored in fastq files, that means that QIIME 2 needs to keep track of it alongside the data. That's where our QIIME 2 artifacts come in. Remember, these are just `.zip` files with a different extension (`.qza`). They can store the fastq data (in the `data/` directory), but also metadata that explicitly defines how quality scores are encoded. When a user imports fastq files into QIIME 2,they must tell QIIME 2 what encoding scheme is used and QIIME 2 can then keep track of it from there. Importing data is when you provide this information to QIIME 2.

## What can be imported in QIIME 2? 



```
qiime tools import
```

Importable types

Importable formats