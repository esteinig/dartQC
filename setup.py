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
      packages=["dartqc"],
      scripts=['bin/dartqc'],
      package_dir={"dartqc": "dartqc"},
      package_data={"dartqc": ["env/dartqc.yaml"]})
