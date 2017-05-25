import os
import json
import pandas
import operator

from dartqc.DartUtils import stamp


class DartPreparator:

    """
    Class for guessing specifications of input file and converting to input for DartReader

    limit: number of beginning rows to read that likely contain the meta data, usually no more than 30

    """

    def __init__(self, file_path, output_name="dartqc", output_path=os.getcwd(), excel_sheet=""):

        self.config = {
            "clone_column": [
                "CloneID"
            ],
            "allele_column": [
                "AlleleID"
            ],
            "sequence_column": [
                "AlleleSequence"
            ],
            "replication_column": [
                "RepAvg"
            ]
        }

        self.file_path = file_path
        self.output_name = output_name
        self.output_path = output_path

        self.excel_sheet = excel_sheet

        self.limit = 30

        self.top = None

        self.data_row = 0
        self.sample_row = 0
        self.data_column = 0

        self.header = None
        self.scheme = dict()

        if self.excel_sheet:
            self._convert_excel()

        self._read_csv()
        self._get_row_indices()
        self._get_column_indices()
        self._reindex()
        self._write_scheme()

    def _convert_excel(self):

        stamp("Converting from Excel")
        stamp("File is", self.file_path)
        stamp("Sheet is", self.excel_sheet)
        data_xls = pandas.read_excel(self.file_path, self.excel_sheet, index_col=None)
        name, ext = os.path.splitext(os.path.basename(self.file_path))
        outfile = os.path.join(self.output_path, name + ".csv")

        stamp("Writing to file", outfile)
        data_xls.to_csv(outfile, encoding='utf-8', index=False)

        self.file_path = outfile

    def _read_csv(self):

        stamp("Loading file", self.file_path)

        self.top = pandas.read_csv(self.file_path, header=None, nrows=30)

    def _get_row_indices(self):

        """

        This function gets the row indices for samples and data, as well as the header, assuming:

            - header row begins after rows starting with "*"
            - data row starts after header row
            - samples are specified in the header row (above calls or raw counts)

        """

        stamp("Guessing data configuration:")

        for i, row in self.top.iterrows():
            if row[0] != "*":
                self.header = row
                self.sample_row = i
                self.data_row = i+1
                self.scheme["sample_row"] = self.sample_row
                self.scheme["data_row"] = self.data_row
                break

    def _get_column_indices(self):

        """

        Get column indices for key words defined in the config file (excel_scheme.json), also get
        column indices for start of data assuming:

            - data starts in first column after placeholder ("*")
            - guessed from the row above the header (self.sample_row-1)

        """

        i = 0
        for col in self.top:
            data = self.top[col]
            if data[self.sample_row-1] != "*":
                self.scheme["data_column"] = i
                break
            i += 1

        for column, aliases in self.config.items():
            column_indices = [pandas.Index(self.header).get_loc(header) for header in aliases]

            if len(column_indices) < 1:
                raise ValueError("Could not find column header", column, "in dataset header. "
                                 "Please provide a different alias in dartQC/schemes/excel_scheme.json.")
            elif len(column_indices) >= 2:
                raise ValueError("Found more than one column header out of:", aliases, "make sure column "
                                 "headers are not repeated in data. You can change aliases that this command is "
                                 "looking for in dartQC/schemes/excel_scheme.json")
            else:
                self.scheme[column] = column_indices[0]

    def _reindex(self):

        """ Reindex values for non-pythonic input to DartReader (better for users) """

        self.scheme = {k: int(v+1) for k, v in self.scheme.items()}

        for k, v in sorted(self.scheme.items(), key=operator.itemgetter(1)):
            stamp(k, "=", v)

        stamp("Please check these values in your data to ensure correct input for DartQC.")

    def _write_scheme(self):

        if self.output_name is None:
            name, ext = os.path.splitext(os.path.basename(self.file_path))
            file_name = name + "_scheme.json"
        else:
            file_name = self.output_name + "_scheme.json"

        out_file = os.path.join(self.output_path, file_name)

        stamp("Writing scheme to:", out_file)

        with open(out_file, "w") as outfile:
            json.dump(self.scheme, outfile, indent=4)

