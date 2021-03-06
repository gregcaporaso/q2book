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
  display_name: calysto_bash
  language: calysto_bash
  name: calysto_bash
---
(metadata)=
# Metadata in QIIME 2

Metadata provides the key to gaining biological insight from your data. In QIIME 2, **sample metadata** may include technical details, such as the DNA barcodes that were used for each sample in a multiplexed sequencing run, or descriptions of the samples, such as which subject, time point, and body site each sample came from in a human microbiome time series. **Feature metadata** is often a feature annotation, such as the taxonomy assigned to an amplicon sequence variant (ASV). Sample and feature metadata are used by many plugins, and examples are provided in this and other chapters illustrating how to work with metadata in QIIME 2.

Metadata is usually specific to a given microbiome study, and compiling sample metadata is typically a step you will have started before beginning your QIIME 2 analysis. It is up to the investigator to decide what information is collected and tracked as metadata. QIIME 2 does not place restrictions on what types of metadata are expected to be present; there are no enforced "metadata standards". This is your opportunity to track whatever information you think may be important to your analyses. When in doubt, collect as much metadata as possible, as you may not be able to retroactively collect certain types of information.

While QIIME 2 does not enforce standards for what types of metadata to collect, the MIxS and MIMARKS standards {cite}`Yilmaz2011-ra` provide recommendations for microbiome studies and may be helpful in determining what information to collect in your study. If you plan to deposit your data in a data archive (e.g. [ENA](https://www.ebi.ac.uk/ena) or [Qiita](https://qiita.ucsd.edu/)), it is also important to determine the types of metadata expected by that resource. Different data archives have their own requirements.

This chapter presents formatting requirements for QIIME 2 metadata files, discussion of how to validate your metadata files for use with QIIME 2, and examples of how metadata can be viewed and used in QIIME 2. 

````{margin}
```{admonition} Attribution
The text in this chapter was imported to this book on 16 March 2021 from content originally written for the QIIME 2 online documentation. The history of the authorship of and other contributions to the original version of this content can be reviewed [here](https://github.com/qiime2/docs/commits/master/source/tutorials/metadata.rst).
```
````

````{margin}
```{admonition} Video
[This video](https://www.youtube.com/watch?v=hh6pqmzJWds) on the QIIME 2 YouTube channel presents a discussion of sample metadata.
```
````

```{admonition} Jargon: metadata files or mapping files? 
You may sometimes hear TSV metadata files referred to as **mapping files**. In the QIIME 1 documentation we often referred to metadata files as mapping files, but we now prefer the term metadata files because it's more descriptive. In QIIME 2, we always refer to these files as **metadata files**, but they are conceptually the same thing as QIIME 1 mapping files. QIIME 2 metadata files are backwards-compatible with QIIME 1 mapping files, meaning that you can use existing QIIME 1 mapping files in QIIME 2 without needing to make modifications to the file.
```

## Metadata Formatting Requirements

QIIME 2 metadata is most commonly stored in a [TSV (i.e. tab-separated values)](https://en.wikipedia.org/wiki/Tab-separated_values) file. These files typically have a ``.tsv`` or ``.txt`` file extension, though it doesn't matter to QIIME 2 what file extension is used. TSV files are simple text files used to store tabular data, and the format is supported by many types of software. TSV files can be edited, imported to, and exported from most spreadsheet programs and databases. Thus, it's usually straightforward to manipulate QIIME 2 metadata using the software of your choosing. If in doubt, we recommend using a spreadsheet program such as Google Sheets to edit and export your metadata files.

Because metadata files contain tabular data, we describe their formatting in terms of **rows** and **columns**. The commonality across QIIME 2 metadata files is that the first [non-comment, non-empty](comments-and-empty-rows) row of the file defines the column headers, and the first column contains a unique identifier for each metadata entry. The following sections describe the formatting requirements for QIIME 2 metadata files. Since there is no universal standard for TSV files, it is important to adhere to these requirements and understand how QIIME 2 will interpret the file's contents to get the most out of your metadata!

```{warning}
Spreadsheet editors often have auto-correct or auto-format features that will modify your data without alerting you that changes will be made {cite}`Ziemann2016-tk`. This is something that you need to watch out for when working with your metadata files in spreadsheet editors.
```

```{tip}
In addition to TSV files, some QIIME 2 Artifacts (i.e. ``.qza`` files) can also be used as metadata. See the section {ref}`artifacts-as-metadata` below for details.
```

(comments-and-empty-rows)=
### Comments and Empty Rows

Rows whose first cell begins with the pound sign (``#``) are interpreted as comments and may appear anywhere in the file. Comment rows are ignored by QIIME 2 and are for informational purposes only. Inline comments (i.e., comments that begin part-way through a row or at the end of a row) are not supported.

Empty rows (e.g. blank lines or rows consisting solely of empty cells) may appear anywhere in the file and are ignored.

(identifier-column)=
### Identifier Column

The first column in the metadata file is the **identifier (ID) column**. This column defines the sample or feature IDs associated with your study. It is not recommended to mix sample and feature IDs in a single metadata file; keep sample and feature metadata stored in separate files.

The **ID column name** (also referred to as the ID column header) must be one of the following values. The values listed below are reserved for use as ID column names and may not be used as IDs or names of other columns in the metadata file.

Case-insensitive (i.e., uppercase or lowercase, or a mixing of the two, is allowed):

- ``id``
- ``sampleid``
- ``sample id``
- ``sample-id``
- ``featureid``
- ``feature id``
- ``feature-id``

````{margin}
```{note}
The case-sensitive ID headers are available for backwards-compatibility with QIIME 1, biom-format, and Qiita files.
```
````

Case-sensitive (i.e., these must appear exactly as presented here):

- ``#SampleID``
- ``#Sample ID``
- ``#OTUID``
- ``#OTU ID``
- ``sample_name``

The following rules apply to IDs:

- IDs may consist of any Unicode characters, with the exception that IDs must not start with the pound sign (``#``), as those rows would be interpreted as comments and ignored. See the section {ref}`identifier-recommendations` for recommendations on choosing identifiers in your study.
- IDs cannot be empty (i.e. they must consist of at least one character).
- IDs must be unique (exact string matching is performed to detect duplicates).
- At least one ID must be present in the file.
- IDs cannot be any of the reserved ID headers listed above.

### Metadata Columns

The ID column is the first column in the metadata file, and can optionally be followed by additional columns defining metadata associated with each sample or feature ID. Metadata files are not required to have additional metadata columns, so a file containing only an ID column is a valid QIIME 2 metadata file. 

The following rules apply to column names:

- May consist of any Unicode characters.
- Cannot be empty (i.e., column names must consist of at least one character).
- Must be unique (exact string matching is performed to detect duplicates).
- Column names cannot use any of the reserved ID headers described in the section {ref}`identifier-column`.

The metadata file line containing the ID column name and any other column names is referred to as the **header row**.

```{admonition} Jargon: metadata columns or metadata categories? 
In previous versions of QIIME 2 and in QIIME 1, *metadata columns* were often referred to as *metadata categories*. Now that we support metadata column typing, which allows you to say whether a column contains *numeric* or *categorical* data, we would end up using terms like *categorical metadata category* or *numeric metadata category*, which can be confusing. We now avoid using the term *category* unless it is used in the context of *categorical* metadata. We've done our best to update our software and documentation to use the term *metadata column* instead of *metadata category*, but there may still be lingering usage of the previous terms out there.
```

### Metadata values

The contents of a metadata file following the ID column and header row (excluding comments and empty lines) are referred to as the **metadata values**. A single metadata value, defined by an (ID, column) pair, is referred to as a **cell**. 

The following rules apply to metadata values and cells:

- May consist of any Unicode characters.
- Empty cells represent *missing data*. Other values such as ``NA`` are not interpreted as missing data; only the empty cell is recognized as "missing". Note that cells consisting solely of whitespace characters are also interpreted as *missing data* because [leading and trailing whitespace characters are always ignored](metadata:whitespace), effectively making the cell empty.

```{note}
The empty cell simply indicates that data is missing, but doesn't indicate what type of missing data it might be. You can use other values of your choosing to denote different types of missing data (e.g. "not applicable" vs. "not collected"). These custom values won't be interpreted as missing data in QIIME 2, but you can still record and use these "missing" metadata values to perform filtering on your data prior to further analyses (e.g. using ``qiime feature-table filter-samples`` to filter samples based on custom "missing" values).
```

(metadata:whitespace)=
### Leading and trailing whitespace characters

If **any** cell in the metadata contains leading or trailing whitespace characters (e.g. spaces, tabs), those characters will be ignored when the file is loaded. Thus, leading and trailing whitespace characters are not significant, so cells containing the values ``'gut'`` and ``'  gut  '`` are equivalent. This rule is applied before any other rules described in this section.

(identifier-recommendations)=
### Recommendations for Identifiers

Our goal with QIIME 2 is to support arbitrary Unicode characters in all cells of metadata files. However, given that QIIME 2 plugins and interfaces can be developed by anyone, we can't make a guarantee that arbitrary Unicode characters will work with all plugins and interfaces. We can therefore make recommendations to users about characters that should be safe to use in identifiers, and we are preparing resources for plugin and interface developers to help them make their software as robust as possible. As developer resources become available, we will announce them in the [Developer Discussion category](https://forum.qiime2.org/c/dev-discussion) on the QIIME 2 Forum.

Sample and feature identifiers with problematic characters tend to cause the most issues for our users. Based on our experiences with QIIME 1, QIIME 2, and other bioinformatics and command line tools, we recommend the following attributes for identifiers:

- Identifiers should be 36 characters long or less.
- Identifiers should contain only ASCII alphanumeric characters (i.e. in the range of ``[a-z]``, ``[A-Z]``, or ``[0-9]``), the period (``.``) character, or the dash (``-``) character.

````{margin}
```{note}
The length recommended here (36 characters or less) is designed to be as short as possible while still supporting version 4 UUIDs formatted with dashes.
```
````

An important point to remember is that sometimes values in your sample metadata can become identifiers. For example, taxonomy annotations can become feature identifiers following ``qiime taxa collapse``, and sample or feature metadata values can become identifiers after applying ``qiime feature-table group``. If you plan to apply these or similar methods where metadata values can become identifiers, you will be less likely to encounter problems if the values adhere to these identifier recommendations as well.

```{tip}
We recommend the [cual-id](https://github.com/johnchase/cual-id) software for assistance with creating sample identifiers. The cual-id paper {cite}`Chase2016-wl` also provides some discussion on how to design identifiers.
```

```{note}
Some bioinformatics tools may have more restrictive requirements on identifiers than the recommendations that are outlined here. For example, Illumina sample sheet identifiers cannot have `.` characters, while we do include those in our set of recommended characters. Similarly, [Phylip](http://evolution.genetics.washington.edu/phylip.html) requires that identifiers are a maximum of 10 characters, while we recommend length 36 or less. If you plan to export your data for use with other tools that may have more restrictive requirements on identifiers, we recommend that you adhere to those requirements in your QIIME 2 analyses as well, to simplify subsequent processing steps.
```

### Column Types

QIIME 2 currently supports *categorical* and *numeric* metadata columns. By default, QIIME 2 will attempt to infer the type of each metadata column: if the column consists only of numbers or missing data, the column is inferred to be *numeric*. Otherwise, if the column contains any non-numeric values, the column is inferred to be *categorical*. Missing data (i.e. empty cells) are supported in categorical columns as well as numeric columns.

QIIME 2 supports an **optional comment directive** to allow users to explicitly state a column's type. This bypasses the column type inference described above. This can be useful if there is a column that appears to be numeric, but should actually be treated as categorical metadata (e.g. a ``Subject`` column where subjects are labeled ``1``, ``2``, ``3``). Explicitly declaring a column's type also makes your metadata file more descriptive because the intended column type is included with the metadata, instead of relying on software to infer the type (which isn't always transparent).

You can use an optional *comment directive* to declare column types in your metadata file. The comment directive must appear *directly* below the header row. The value in the ID column in this row must be ``#q2:types`` to indicate the row is a comment directive. Subsequent cells in this row may contain the values ``categorical`` or ``numeric`` (both case-insensitive). The empty cell is also supported if you do not wish to assign a type to a column (the type will be inferred in that case). Thus, it is easy to include this comment directive without having to declare types for every column in your metadata.

````{margin}
```{tip}
The command ``qiime metadata tabulate`` can be used to review the column types of your QIIME 2 Metadata. This works whether you're using the comment directive, type inference, or a combination of the two approaches.
```
````

````{margin}
```{note}
The ``#q2:types`` comment directive is the only supported comment directive; others may be added in the future (e.g. ``#q2:units``). For this reason, rows starting with ``#q2:`` are disallowed, as we reserve that namespace for future comment directives.
```
````

### Number Formatting

If a column is to be interpreted as a *numeric* metadata column (either through column type inference or by using the ``#q2:types`` comment directive), numbers in the column must be formatted following these rules:

- Use the decimal number system: ASCII characters ``[0-9]``, ``.`` for an optional decimal point, and ``+`` and ``-`` for positive and negative signs, respectively.

  - Examples: ``123``, ``123.45``, ``0123.40``, ``-0.000123``, ``+1.23``

- Scientific notation may be used with *E-notation*; both ``e`` and ``E`` are supported.

  - Examples: ``1e9``, ``1.23E-4``, ``-1.2e-08``, ``+4.5E+6``

- Only up to 15 digits **total** (including before and after the decimal point) are supported to stay within the 64-bit floating point specification. Numbers exceeding 15 total digits are unsupported and will result in undefined behavior.

- Common representations of *not a number* (e.g. ``NaN``, ``nan``) or infinity (e.g. ``inf``, ``-Infinity``) are **not supported**. Use an empty cell for missing data (e.g. instead of ``NaN``). Infinity is not supported at this time in QIIME 2 metadata files.

### Advanced File Format Details

If you're creating TSV files by manually (e.g. in a text editor, not a spreadsheet program such as Microsoft Excel or Google Sheets) or writing your own software to consume or produce QIIME 2 metadata files refer to [](advanced-metadata-formatting). 

(metadata:validation)=
## Metadata Validation

QIIME 2 will automatically validate a metadata file anytime it is used. This will inform you of any errors in your metadata formatting, which you can then correct. 

Sample and feature metadata files stored in Google Sheets can additionally be validated using [Keemei](https://keemei.qiime2.org/). Keemei can be installed as a Google Sheets Add-on following the instructions on the Keemei website. After installing the Keemei Add-on, select *Add-ons > Keemei > Validate QIIME 2 metadata file* to validate metadata stored in Google Sheets.

Using [Keemei](https://keemei.qiime2.org/) to validate your metadata is typically more convenient than using QIIME 2's built-in validation because a report of all validation errors and warnings will be presented each time Keemei is run. [Keemei](https://keemei.qiime2.org/) also will warn users about identifiers that don't meet the [ID recommendations](identifier-recommendations) presented above. Loading your metadata in QIIME 2 will typically present only a single error at a time, which can make identifying and resolving validation issues cumbersome, especially if there are many issues with the metadata. However, if you're unable to use Keemei (for example because of firewall issues that disallow accessing Google services, or due to concerns over security), QIIME 2's built-in validation is totally sufficient. 

````{margin}
```{note}
We ultimately plan to support Keemei-style validation within QIIME 2.
```
````

## Using Metadata Files

In this section you'll learn how to use metadata in QIIME 2 by looking at some examples. To follow along start by downloading an example sample metadata TSV file:

```{code-cell}
:tags: [hide-output]
curl -sL \
  "https://data.qiime2.org/2021.4/tutorials/moving-pictures/sample_metadata.tsv" > \
  "sample-metadata.tsv"
```

Since this is a TSV file, it can be opened and edited in a variety of applications, including text editors, Microsoft Excel, and Google Sheets (e.g. if you plan to validate your metadata with [Keemei](https://keemei.qiime2.org/)).

QIIME 2 also provides a visualizer for viewing metadata in an interactive table:

```{code-cell}
qiime metadata tabulate \
    --m-input-file sample-metadata.tsv \
    --o-visualization tabulated-sample-metadata.qzv
```

Load the resulting file, as well as other QIIME 2 visualizations generated in this chapter, using [QIIME 2 View](https://view.qiime2.org).

````{margin}
```{admonition} Video
[This video](https://t.co/eJbm03cnSa) on the QIIME 2 YouTube channel illustrates how to use QIIME 2 View. 
```
````

```{admonition} Exercise
Based on the table in ``tabulated-sample-metadata.qzv``, how many samples are associated with ``subject-1``? How many samples are associated with the ``gut`` body site? Hint: use the search box and/or the column sorting options to assist with this query.
```

(artifacts-as-metadata)=
## Using QIIME 2 Artifacts as Metadata

In addition to TSV metadata files, QIIME 2 also supports viewing some kinds of artifacts as metadata. An example of this is artifacts of type ``SampleData[AlphaDiversity]``.

To get started with understanding artifacts as metadata, first download an example artifact:

```{code-cell}
curl -sL \
  "https://data.qiime2.org/2021.4/tutorials/metadata/faith_pd_vector.qza" > \
  "faith_pd_vector.qza"
```

To view this artifact as metadata, simply pass it in to any method or visualizer that expects to see metadata (e.g. ``metadata tabulate`` or ``emperor plot``):

```{code-cell}
qiime metadata tabulate \
    --m-input-file faith_pd_vector.qza \
    --o-visualization tabulated-faith-pd-metadata.qzv
```

```{admonition} Exercise
What is the largest value of Faith's PD? What is the smallest? Hint: use the column sorting functions to assist with this query.
```

When an artifact is viewed as metadata, the result includes that artifact's provenance in addition to its own.

```{admonition} Exercise
Identify this artifact in the interactive provenance graph after loading the `tabulated-faith-pd-metadata.qzv` file with [QIIME 2 View](https://view.qiime2.org).
```

## Merging metadata

Since metadata can come from many different sources, QIIME 2 supports metadata merging when running commands. Building upon the examples above, simply passing ``--m-input-file`` multiple times will combine the metadata columns in the specified files:

```{code-cell}
qiime metadata tabulate \
    --m-input-file sample-metadata.tsv \
    --m-input-file faith_pd_vector.qza \
    --o-visualization tabulated-combined-metadata.qzv
```

The resulting metadata after the merge will contain the intersection of the identifiers across all of the specified files. In other words, the merged metadata will only contain identifiers that are shared across all provided metadata files. This is an *inner join* using database terminology.

```{admonition} Exercise
Modify the command above to merge the [evenness vector](https://docs.qiime2.org/2021.4/data/tutorials/moving-pictures/core-metrics-results/evenness_vector.qza) of ``SampleData[AlphaDiversity]`` after the Faith's PD vector. What happens when merging the three artifacts? How many columns are present in the resulting metadata visualization? How many of those columns represent the sample IDs? How many of those columns represent ``SampleData[AlphaDiversity]`` metrics? What happens to the visualization if the order of the metadata files is reversed? Hint, take a closer look at the column ordering.
```

Metadata merging is supported anywhere that metadata is accepted in QIIME 2. For example, it might be interesting to color an Emperor plot based on the study metadata, or sample alpha diversity. This can be accomplished by providing both the sample metadata file *and* the ``SampleData[AlphaDiversity]`` artifact:

```{code-cell}
curl -sL \
  "https://data.qiime2.org/2021.4/tutorials/metadata/unweighted_unifrac_pcoa_results.qza" > \
  "unweighted_unifrac_pcoa_results.qza"

qiime emperor plot \
    --i-pcoa unweighted_unifrac_pcoa_results.qza \
    --m-metadata-file sample-metadata.tsv \
    --m-metadata-file faith_pd_vector.qza \
    --o-visualization unweighted-unifrac-emperor-with-alpha.qzv
```

```{admonition} Exercise
What body sites are associated with the highest Faith's phylogenetic diversity value? Hint: first color by body site, and then color by Faith's PD using a continuous color scheme.
```


(exploring feature metadata)=
## Exploring feature metadata

Metadata in QIIME 2 can be applied to sample or features --- so far we have only dealt with sample metadata. This section will focus on feature metadata, specifically how to view ``FeatureData`` as metadata.

To get started with feature metadata, first download the example files:

```{code-cell}
curl -sL \
  "https://data.qiime2.org/2021.4/tutorials/metadata/rep-seqs.qza" > \
  "rep-seqs.qza"

curl -sL \
  "https://data.qiime2.org/2021.4/tutorials/metadata/taxonomy.qza" > \
  "taxonomy.qza"
```

We have downloaded a ``FeatureData[Sequence]`` file (``rep-seqs.qza``) and a ``FeatureData[Taxonomy]`` file (``taxonomy.qza``). We can merge (and ``tabulate``) these files to associate the representative sequences with their taxonomic annotations:

```{code-cell}
qiime metadata tabulate \
    --m-input-file rep-seqs.qza \
    --m-input-file taxonomy.qza \
    --o-visualization tabulated-feature-metadata.qzv
```

The resulting table shows the joined metadata files with a column of the the feature IDs, a column of the representative sequences, a column of the taxonomic assignments, and lastly, a column of the assignment confidence.

```{admonition} Exercise
Are all artifacts (``.qza`` files) viewable as metadata? Hint: try tabulating a [`FeatureTable` artifact](https://docs.qiime2.org/2021.4/data/tutorials/moving-pictures/table.qza). Are all metadata files stored as ``.qza`` files?
```

Finally, there are export options available in the visualizations produced from ``metadata tabulate``. Using the results from ``tabulated-feature-metadata.qzv``, export the data as a new TSV. Open that file in a TSV viewer or text editor and note that the contents are the same as the interactive metadata table in the visualization.

```{admonition} Exercise
Can the exported TSV from the above step be used as metadata? What are some benefits of being able to export metadata (hint: see the discussion above about metadata merging)? What about some potential drawbacks (hint: what happens to data provenance when data is exported from QIIME 2)?
```

## List of works cited

```{bibliography} ../references.bib
:filter: docname in docnames
```
