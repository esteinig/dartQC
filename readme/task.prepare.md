## Task: Prepare

This task's main function is to generate a scheme file, so that subsequent modules know where to find the right rows and columns in `--file`.  **Input is the double-row format for SNPs by DArT**. Executing this task will attempt to guess which rows and columns the data are in. You can also specify a sheet name with `--sheet`, which will convert an Excel sheet from `--file` to CSV. This task needs to beb run for both raw and called read files, if you are later using task `process`.

If file (CSV) is in current working directory:

`dartqc.py prepare --file example.csv`

This produces output: `example_scheme.json`

If file (Excel) with spreadsheet name (double_row) in current working directory:

`dartqc.py prepare --file example.xlsc --sheet double_row`

This produces outputs: `example_scheme.json`, `example.csv`

You can change the output directory with the global option `-o` and change the scheme file name with the task option `--name`:

`dartqc.py -o ./prep prepare --file example.csv --name prep`

This produces output: `./prep/prep_scheme.json`

### Data Scheme: Assumptions

### Data Scheme: Manual

