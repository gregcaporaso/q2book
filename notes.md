# Notes to self...

## Building the book and publishing to GH pages

From the top-level directory of the repository:

```
$ jupyter-book build book
```

Following the instructions presented in the Jupyter Book documentation [here](https://jupyterbook.org/publish/gh-pages.html), I publish to GH Pages after building as follows:

```
$ ghp-import -n -p -f book/_build/html
```
