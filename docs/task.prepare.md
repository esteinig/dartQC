## Task: Prepare

This task's main function is to generate a scheme file, so that subsequent modules know where to find the right rows and columns in the input data `--file`.  **Input is the double-row format for SNPs by DArT**. 

Executing this task will attempt to guess which rows and columns the data are in and output a JSON. You can also specify a sheet name with `--sheet`, which will convert an Excel sheet from `--file` to CSV. This task needs to be run for both raw and called read files to generate both `--raw_scheme` and `--call_scheme`, if you are later using the task [`process`](https://github.com/esteinig/dartQC/blob/master/readme/task.process.md).

---

If file (CSV) is in current working directory:

`dartqc.py prepare --file example.csv`

This produces output: `example_scheme.json`

---

If file (Excel) with spreadsheet name (double_row) in current working directory:

`dartqc.py prepare --file example.xlsc --sheet double_row`

This produces outputs: `example_scheme.json`, `example.csv`

---

You can change the output directory with the global option `-o` and change the scheme file name with the task option `--name`:

`dartqc.py -o ./prep prepare --file example.csv --name prep`

This produces output: `./prep/prep_scheme.json`

---

### Data Formatting

Note that the data from DArT needs to include `CloneID` as a simple number, such as: `12753571`.

### Data Scheme: Manual

You can also reate your own manual file if the task fails to prepare the right format for you, simply create a JSON file like this and supply the correct rows and columns (without the comments):

```
{
  "clone_column": 1,                # CloneID
  "allele_column": 2,               # AlleleID
  "sequence_column": 3,             # AlleleSequence
  "replication_column": 17,         # RepAvg
  "data_column": 18,                # Column, start of read or call data
  "sample_row": 7,                  # Row, contains sample designations
  "data_row": 8                     # Row, start of data (SNPs)
}
```

This file can then be specified in tasks [`process`](https://github.com/esteinig/dartQC/blob/master/readme/task.process.md) and [`filter`](https://github.com/esteinig/dartQC/blob/master/readme/task.filter.md) via `--raw_scheme` or `--call_scheme`.

### Data Scheme: Assumptions

The task scans the first thirty rows of the input file for the following:

- header row starting in the first row that does not start with `*`
- sample ids in the same row as header, i.e. above the data
- data starts one row after the header row
- looks for columns "AlleleSequence", "CloneID", "AlleleID" and "RepAvg" in header
- data and sample ids start in first column after `*` in the row above the header
- output are non-pythonic indices (starting with 1)

See example files in this repository for the correct format.
