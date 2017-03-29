import textwrap


class DartMessages:

    def __init__(self):

        pass

    def get_filter_message(self, filter, threshold, initial, removed, retained):

        filter_msg = textwrap.dedent("""
                SNP Filter
        -------------------------------

        {0} at {1}

        Initial:    {2}
        Removed:    {3}
        Retained:   {4}

        -------------------------------
        """ .format(filter.upper(), threshold, initial, removed, retained))

        return filter_msg

    def get_redundancy_message(self, type, initial, removed, retained):

        redundancy_msg = textwrap.dedent("""
                  REDUNDANCY
        -------------------------------

        {0}

        Initial:    {1}
        Removed:    {2}
        Retained:   {3}

        -------------------------------
        """ .format(type.upper(), initial, removed, retained))

        return redundancy_msg

    def get_cdhit_message(self, identity):

        cluster_msg = textwrap.dedent("""
                  CLUSTERING
        -------------------------------

        Running CDHIT-EST...

        Threshold: {0}%

        -------------------------------
        """ .format(identity*100))

        return cluster_msg