from dartqc.DartReader import DartReader
from dartqc.DartModules import SampleModule, PopulationModule
import os

dart_reader = DartReader()
dart_reader.set_options(project="test", out_path=os.getcwd(), scheme="prawn_data_scheme.json")

dart_reader.read_double_row(file="prawn_data.csv", basic=True)

dart_reader.read_pops("prawn_pops.csv", sep=",")

data, attributes = dart_reader.get_data()

pm = PopulationModule(data, attributes)
data, attributes = pm.get_data(mono="all")

#im = SampleModule(data, attributes)

#im.filter_data()