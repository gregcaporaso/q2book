(demultiplexing)=
# Demultiplexing a sequencing run

**TODO: All content is currently copy/pasted from its previous home in the overview tutorial chapter. I decided to simplify that chapter as much as possible. This content will need to be contextualized for this chapter.**

The second file that we downloaded above was our multiplexed sequence data. This is sequencing data as it came off an Illumina DNA sequencing instrument. At this stage, the sequences are all grouped together - they are not grouped by sample. Instead, during sample preparation (PCR specifically), a short DNA sequence referred to as an index or a barcode was added to sequences on a per-sample basis. For example, `sample-1` might be assigned the barcode `TGACCGTACGTA`, and `sample-2` might be assigned the barcode `TGGTAGACCCGT`. All sequences derived from `sample-1` will have the `sample-1` barcode associated with them, and all sequences derived from `sample-2` will have the `sample-2` barcode associated with them. One piece of sample metadata in our sample metadata file is the barcode sequence associated with each sample. Assign sequences to the sample they are derived from is referred to as _demultiplexing_ the sequencing run. The `q2-demux` plugin in QIIME 2 has several methods for demultiplexing sequencing runs. These take the sequencing data as input, as well as the sample metadata. 

````{margin}
```{admonition} Jargon
_Demux_ is a common abbreviation for _demultiplex_. For example, if someone asks if your sequence data is demuxed, they're probably wondering whether the sequences from your samples are all group together (i.e., they are multiplexed), or if they have already been associated with the samples they are derived from (i.e., they are demultiplexed).
```
````

````{margin}
```{note}
You can start a QIIME 2 analysis with data is that is multiplexed, as illustrated here, or data that is demultiplexed. Personally I prefer to start with data that is already demultiplexed - I typically ask the group doing the sequencing for me to demultiplex the data before delivering it. Because different sequencing centers use different multiplexing strategies (e.g., single barcoding versus dual barcoding) I find that this is easier because I don't have to figure out how they barcoded the data, and then figure out how to demultiplex it. Since they know how they barcoded the data, they'll know how to demultiplex it. This is still a bit of a personal preference though: one option is not necessarily better than the other, and you can start your analysis with whatever you have. You'll probably develop your own preference as you gain experience.
```
````

## Demultiplexing ...

As a first analysis step with QIIME 2, demultiplex the sequence data with the following command:

```
qiime demux emp-paired \
  --i-seqs emp-paired-end-sequences.qza \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --p-rev-comp-mapping-barcodes \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details demux-log.qza
```

Since this is one of the first QIIME 2 commands that you've run (at least while reading this book) let's take a look at it in some detail. As we look at this command, I'm going to assume that you either have some command line software experience, or that you've read _the introductory chapter in this book on using the command line_ (**TODO: this chapter doesn't exist yet.**). There are four components to notice on the first line of this command, and as is typically the case when running command line software, the command's components are separated by spaces. The first space-separated component is `qiime`. This tells your computer that the QIIME 2 program should be run. The second space-separated component is `demux`. This is the name of the q2-demux plugin, as far as QIIME 2 is concerned, so this tells QIIME 2 that you want to use the `q2-demux` plugin. The third space-separated component is `emp-paired`, which is an action in the q2-demux plugin. This action demultiplexes paired-end sequencing data that is formatted as is typical for Earth Microbiome Project (EMP) sequenced data (more on that later). As mentioned in the previous chapter, if you'd like to learn about this action, you could run the command `qiime demux emp-paired --help`, which will print help text to the screen. Finally, the `\` at the end of the line is worth mentioning now. When you're working on the command line, line breaks (i.e., the character that is received when you press the _Return_ key on your keyboard) signifies that you have finished entering the command and that the command should now be executed (i.e., run). Since QIIME 2 commands can be long, it's helpful when documenting them to split them over multiple lines so that they can be read without the reader having to scroll to the right. If the documentation is being presented in a non-interactive medium, such as a printed book, the formatting of the command could actually be misleading - for example, it could look like the command should be split across multiple lines. Splitting long commands across multiple lines also therefore helps to ensure accuracy of the documentation and it improves readability. The `\` character here simply means that the terminal shouldn't interpret the line break that follows it as the end of the command, but rather that the command will continue on the next line. The same command, without the line breaks, would look like the following:

``` 
qiime demux emp-paired --i-seqs emp-paired-end-sequences.qza --m-barcodes-file sample-metadata.tsv --m-barcodes-column barcode-sequence --p-rev-comp-mapping-barcodes --o-per-sample-sequences demux.qza --o-error-correction-details demux-log.qza
```

That command will do the exact same thing as the one above. I think you'll agree though that the first variant, where the command is split across multiple lines, is easier to read. 

After we specify the program (`qiime`), the plugin (`demux`), and the action (`emp-paired`) that we want to run, we provide options to that action to tell it what files we want it to operate on, what we want it to do exactly, and what we want the outputs to be called. These components of the command are specified on the remaining lines. Let's look at those now. 

This action takes one QIIME 2 artifact as input, the multiplexed sequences that we downloaded above. In the QIIME 2 command line interface, input artifacts are always specified with options beginning with `--i-`. This action also takes a metadata file and a metadata column, both of which are specified with options beginning with `--m-` (where the _m_ implies _metadata_). Here, these specify the file containing the metadata (that's the second file that we downloaded above) and column header in that file that contains the barcode sequences. If you have opened the `sample-metadata.tsv` file, you'll see that there is a column called `barcode-sequence`. This column could be called nearly anything (there are some rules for what column names are allowed that we'll cover later). As long as you specify the column name correctly when you call this action, it will know where to find the barcodes. Input and metadata parameters, in general, are files that you're providing to a QIIME 2 action for actions to be performed on. They will never be modified by QIIME 2 - they are only viewed or read. In this command we're also providing one parameter, specified by the option beginning with `--p-`. Parameters modify the behavior of an action, and parameters that don't take any values (such as this one) are often referred to as _flags_. This particular parameter tells the action that the orientation of the barcodes in the sample metadata file are the reverse complement of the barcodes in the sequence data, and thus that the barcodes from the sample metadata should be reverse complemented before being used for demultiplexing. We'll see examples of parameters that take a value later in this chapter. Finally, this command will generate two files as output, and we specify the paths (absolute or relative) where we want to store those files with output options. Output options always begin with `--o-` in QIIME 2. The first output created by this command was `demux.qza`, which contains our demultiplexed sequences (or sequences that have been grouped by sample). The second output created by this command is `demux-log.qza`. This is a log file that contains some details about this action. 

