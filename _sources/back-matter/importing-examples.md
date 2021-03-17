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
(importing-examples)=
# Appendix: Examples illustrating importing data into QIIME 2

This document will ultimately pull content from the [QIIME 2 documentation's importing guide](https://docs.qiime2.org/2020.11/tutorials/importing/). Refer there for now. 

(importing-examples:available-semantic-types)=
## Discovering the available importable semantic types

Plugins define semantic types, so the semantic types that can be imported depends on what plugins you have installed. To discover what semantic types can be imported in your QIIME 2 deployment, you can run the following command:

```{code-cell}
qiime tools import --show-importable-types
```

## Discovering the available importable formats

Plugins define file types (or formats) as well. To discover what formats can be imported from in your QIIME 2 deployment, you can run the following command:

```{code-cell}
qiime tools import --show-importable-formats
```