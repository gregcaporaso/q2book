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

# Reading this book

This book is intended to be read sequentially and will teach you everything you need to know to get started using QIIME 2 for your own microbiome research using example data sets along the way. As you learn to use the tools, you'll be introduced to underlying theory in microbiome bioinformatics. I also sprinkle in tips and tricks on how to get more out of QIIME 2, including how I use QIIME 2, so you'll be a QIIME 2 power user by the time you're done!  

If you'd like to learn more about the underlying algorithms used in microbiome research, many of those are covered in my other book, _[An Introduction to Applied Bioinformatics 2nd Edition](http://readIAB.org)_. Many of the algorithms discussed in that book are fundamentals of bioinformatics, so if you're interested in becoming a bioinformatician or bioinformatics software developer, you should read these.  

## Reading interactively with Binder

This book can be read statically or interactively. The simplest way to read interactively is with [Binder](https://mybinder.org/). When chapters are available to read interactively, you'll see a rocket ship icon toward the top-right of the page. The current page can be read interactively, allowing you to execute the following code block (making any modifications you choose to make). If you'd like to try it out, click the rocket ship icon, and then the Launch Binder box that pops up beneath it. 

Code blocks will be presented throughout the book. If you're reading a static version of the book, the output of running the code will be presented below the code block. If you're reading an interactive version of the book, you'll have to run the code cell to see the output. 

Here's an example Python 3 code cell that imports the scikit-bio library and prints some information to the screen. 

```{code-cell}
import skbio

print(skbio.title)
print(skbio.art)
```

## Conventions used in this text

### Text formats

I have tried to be consistent in how I emphasize text throughout the book. **Bolded text** indicates new terms that are being introduced; _italicized text_ is used for emphasizing important ideas; and `fixed-width text` is used to indicate literal values, for example names of variables that are being used in the text. 

### Special text blocks

Throughout the text, you'll sometimes see blocks of text that stand out as follows:

```{admonition} Tip: Learn programming a little at a time
:class: tip
Spending a few minutes coding every day is a great way to build your programming skills, empower your own research, and make you more competitive for future career opportunities. Python 3 and R are great languages to get started with for bioinformatics, but just about any language will be fine for learning. 
```

These are intended to share ideas related to what is being presented in the text. I'll use a few types of these special text blocks:

- Video: a video that may help you understand this content. 
- Tip: an idea that may help you in your learning.
- Note: a idea that is related to the content being presented, but is a bit off topic relative to the current discussion.
- Warning: a common error that is encountered or mistake that is made.
- Food for thought: a question or topic to think about related to the current content.
- Exercise: a suggestion for something to experiment with to help you understand the current content.
- Jargon: a term that is common among bioinformatics practitioners, but which may not be immediately clear to beginners. 
- Attribution: used to indicate ideas or content that were derived from other open access resources, such as posts on the QIIME 2 Forum or from the QIIME 2 documentation.
