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

## Reading interactively with Binder

This book can be read statically or interactively. The simplest way to read interactively is with [Binder](https://mybinder.org/). When chapters are available to read interactively, you'll see a rocket ship icon toward the top-right of the page. The current page can be read interactively, allowing you to execute the following code block (making any modifications you choose to make). If you'd like to try it out, click the rocket ship icon, and then the Launch Binder box that pops up beneath it. 

Code blocks will be presented throughout the book. If you're reading a static version of the book, the output of running the code will be presented below the code block. If you're reading an interactive version of the book, you'll have to run the code cell to see the output. 

Here's an example Python 3 code cell that imports the scikit-bio library and prints some information to the screen. 

```{code-cell}
import skbio

print(skbio.title)
print(skbio.art)
```


## Conventions 

Text formats:

- **bolded text** indicates new terms that are being introduced. These will ultimately show up in the glossary as well. 
- _italicized text_ is used for emphasizing import ideas. 
- `fixed-width text` is used to indicate literal values, for example file names or parameters that are used by a command. 

## Special text blocks

Throughout the text, you'll sometimes see blocks of text that stand out as follows:

```{admonition} Tip: Learn programming a little at a time
:class: tip
Spending a few minutes coding every day is a great way to build your programming skills, empower your own research, and make you more competitive for future career opportunities. Python 3 and R are great languages to get started with for bioinformatics, but just about any language will be fine for learning. 
```

These are intended to share ideas related to what is being presented in the text. I'll use a few types of these special text blocks:

- Video: a video is available that may help you understand this content. 
- Tip: an idea that may help you in your learning.
- Warning: a common error that is encountered. 
- Note: a idea that is related to the content being presented, but is a bit off topic relative to the current discussion. 
- Food for thought: a question or topic to think about related to the current content
- Interactive exercise: a suggestion for something to experiment to help you develop your understanding of the current content
- Jargon: discussed a term that is common among microbiome or bioinformatics practitioners, but which may not be immediately clear to beginners. 
