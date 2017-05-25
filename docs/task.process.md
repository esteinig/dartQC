# Task: Process

```
dartqc process [--help] --raw [--raw_scheme] --calls [--call_scheme] [--read_sum]

Arguments:

--raw, -r         path to raw read count csv file
--raw_scheme      path to raw scheme json file
--calls, -c       path to call csv file
--call_scheme     path to call scheme json file
--read_sum        set all calls to missing where sum of read counts < read_sum
```

This tasks runs a pre-processing step on the call data, given raw read counts that can be requested from DArT. At the moment, the pre-processing is based on the sum of both allele counts for each SNP:

1. Inputs are the raw and called data files; the module first checks for congruence of samples and SNP IDs (`AlleleID`), keeping the intersection between raw and called data.
2. Duplicate SNP columns are collapsed in the raw read file and their read counts summed for each allele of the SNP.
3. Both allele's raw read counts are summed for each SNP. If their sum is smaller than `--read_threshold` the call is set to missing (`-`)

Example: SNP with ID `123144124` has 3 total counts for Allele 1 and 4 total counts for Allele 2, their sum is `3 + 4 = 7` and is therefore silenced at default threshold of 10.

Output is the data and its attributes as JSON: `project_data.json` and `project_attr.json`. These files can be passed into task [`filter`](https://github.com/esteinig/dartQC/blob/master/readme/task.filter.md) using the flag `--processed`.

Make sure you have generated the scheme files for both raw and call data manually or with task [`prepare`](https://github.com/esteinig/dartQC/blob/master/readme/task.prepare.md).

---

For call and raw read data in current working directory:

`dartqc.py --project preprocess process --calls example_calls.csv --call_scheme example_calls_scheme.json --raw example_raw.csv --raw_scheme example_raw_scheme.json --read_threshold 10`

This generates outputs: `preprocess_data.json` and `preprocess_attr.json`

---
