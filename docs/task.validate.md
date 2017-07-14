# Task: Validate

```
dartqc validate [-h] --raw RAW_FILE --raw_scheme RAW_SCHEME
                     --calls CALL_FILE --call_scheme CALL_SCHEME
                     --id_list ID_LIST [--cdhit_path CDHIT_PATH]

Arguments:
-h, --help          show this help message and exit
--raw, -r           path to raw read file
--raw_scheme        path to raw scheme json file
--calls, -c         path to called read file
--call_scheme       path to call scheme json file
--id_list, -i       path to CSV file with list of official clone IDs that should be used (eg. to fix ID's that Dart outputs
                    wrong)
--cdhit_path        Path to the cdhit 2d executable (required if cd-hit-2d doesn't work on cmd line)
```

Task to validate the Clone ID's and SNP locations (by clusteringing with cd-hit-est-2d)

Validation information is output to the <project>_seq_vals.csv file and new data and read count files are output with
clone ID's renamed based on the passed in id_list.  These new files are <project>_data_validated.csv and
<project>_read_counts_validated.csv respectively

This task is expected to be run as a stand-alone operation before pre-processing & generates new data and read count files to run the down line processing with.

Note:  If the cdhit_path given is actually for cd-hit-est it will auto rename to cd-hit-est-2d