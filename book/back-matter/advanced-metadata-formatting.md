(advanced-metadata-formatting)=
# Appendix: Advanced metadata formatting

If you're creating TSV files by manually (e.g. in a text editor) or writing your own software to consume or produce QIIME 2 metadata files this section provides additional formatting details. If you're creating and exporting QIIME 2 metadata files using a spreadsheet program (e.g. Microsoft Excel, Google Sheets) you can skip this content. 

## TSV Dialect and Parser

QIIME 2 attempts to interoperate with TSV files exported from Microsoft Excel, as this is the most common TSV "dialect" we have seen in use. The QIIME 2 metadata parser (i.e. reader) uses the [Python csv module](https://docs.python.org/3/library/csv.html) ``excel-tab`` dialect for parsing TSV metadata files. This dialect supports wrapping fields in double quote characters (``"``) to allow for tab, newline, and carriage return characters within a field. To include a literal double quote character in a field, the double quote character must be immediately preceded by another double quote character. See the [Python csv module](https://docs.python.org/3/library/csv.html) for complete documentation on the ``excel-tab`` dialect.

## Encoding and Line Endings

Metadata files must be encoded as UTF-8, which is backwards-compatible with ASCII encoding.

Unix line endings (``\n``), Windows/DOS line endings (``\r\n``), and "classic Mac OS" line endings (``\r``) are all supported by the metadata parser for interoperability. When metadata files are written to disk in QIIME 2, the line endings will always be ``\r\n`` (Windows/DOS line endings).

## Trailing Empty Cells and Jagged Data

The metadata parser ignores any trailing empty cells that occur past the fields declared by the header. This is mainly for interoperability with files exported from some spreadsheet programs. These trailing cells/columns may be jagged (or not); they will be ignored either way when the file is read.

If a row doesn't contain as many fields as declared by the header, empty cells will be padded to match the header length (again, this is mainly for interoperability with exported spreadsheets).
