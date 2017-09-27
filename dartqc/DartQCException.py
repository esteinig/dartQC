import inspect
import sys
import traceback

import re


class DartQCException(RuntimeError):
    def __init__(self, *msg):
        message = ""
        for arg in msg:
            message += str(arg)

        try:
            ln = sys.exc_info()[-1].tb_lineno
            file = traceback.extract_tb(sys.exc_info()[-1])
            if isinstance(file, list):
                file = file[0].filename
        except AttributeError:
            ln = inspect.currentframe().f_back.f_lineno
            file = inspect.getframeinfo(inspect.currentframe().f_back).filename

        if file is not None:
            file = re.compile("[\\\/]+").split(file)
            file = file[-1]

        self.args = "Error ({3}:{1}): {2}".format(type(self), ln, message, file),
        sys.exit(self)

    # def __init__(self, *args, **kwargs):
    #     RuntimeError.__init__(self, args, kwargs)