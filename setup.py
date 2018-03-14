from distutils.core import setup


def readme():
    with open('README.md', 'r') as file:
        file.read()

setup(name='dartqc',
      version='0.1.5',
      description='Quality control for SNP data from Diversity Array Technologies (DArT)',
      url='http://github.com/esteinig/dartqc',
      download_url='http://github.com/esteinig/dartqc/archive/0.1.5.tar.gz',
      author='Eike J. Steinig',
      author_email='eikejoachim.steinig@my.jcu.edu.au',
      license='MIT',
      packages=["dartqc", "dartqc.filters", "dartqc.input", "dartqc.output"],
      scripts=['bin/filter.py', 'bin/install.py', 'bin/CreatePBSScript.py'],
      package_dir={"dartqc": "dartqc", "filters": "dartqc/filters", "input": "dartqc/input", "output": "dartqc/output"},
      package_data={"dartqc": ["env/dartqc.yaml"]}, requires=['simplejson', 'numpy'])
