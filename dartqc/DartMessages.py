import textwrap

from dartqc.DartUtils import stamp

class DartMessages:

    def __init__(self):

        pass

    def get_filter_message(self, filter, threshold, initial, removed, retained, time=True):

        filter_msg = textwrap.dedent("""
                SNP Filter
        -------------------------------

        {0} at {1}

        Initial:    {2}
        Removed:    {3}
        Retained:   {4}

        -------------------------------
        """ .format(filter.upper(), threshold, initial, removed, retained))

        if time:
            stamp("Filtered {0} at {1}".format(filter.upper(), threshold))
            stamp("Removed {0} SNPs".format(removed))
            stamp("Retained {0} SNPs".format(retained))
        else:
            print(filter_msg)

    def get_redundancy_message(self, type, initial, removed, retained, time=True):

        redundancy_msg = textwrap.dedent("""
                  REDUNDANCY
        -------------------------------

        {0}

        Initial:    {1}
        Removed:    {2}
        Retained:   {3}

        -------------------------------
        """ .format(type.upper(), initial, removed, retained))

        if time:
            stamp("Redundancy module {0}".format(type.upper()))
            stamp("Removed {0} SNPs".format(removed))
            stamp("Retained {0} SNPs".format(retained))
        else:
            print(redundancy_msg)

    def get_cdhit_message(self, identity, time=True):

        cluster_msg = textwrap.dedent("""
                  CLUSTERING
        -------------------------------

        Running CDHIT-EST...

        Threshold: {0}%

        -------------------------------
        """ .format(identity*100))

        if time:
            stamp("Running CD-HIT at nucleotide identity {0}%".format(identity*100))
        else:
            print(cluster_msg)