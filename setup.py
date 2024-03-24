from setuptools import setup, find_packages

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

VERSION = '1.0.5' 
DESCRIPTION = 'TnSeeker'
LONG_DESCRIPTION = 'Versatile Python3 module for processing and analyzing anything related with Tn-Seq. Requires bowtie2 to be callable from terminal. File wise, only fastq and genbank annotation files are required.'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="tnseeker", 
        version=VERSION,
        author="Afonso M Bravo",
        author_email="<afonsombravo@hotmail.com>",
        url='https://github.com/afombravo/tnseeker',
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type='text/markdown',
        packages=find_packages(),
        install_requires=["numpy >= 1.19.2",
                           "matplotlib >= 3.3.4",
                           "numba >= 0.53.1",
                           "biopython >= 1.81",
                           "datetime",
                           "argparse",
                           "regex",
                           "multiprocess",
                           "pathlib",
                           "scipy >= 1.6.2",
                           "seaborn >= 0.12.1",
                           "statsmodels >= 0.12.2"],
        
        keywords=['tn-seq', 'essentiality','rb-seq'],
        classifiers= [
            "Development Status :: 3 - Alpha", #change here for when the product is more finished
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3",
            "Operating System :: OS Independent",
        ],
        include_package_data=True,
        package_data={'': ['data/*.csv']},
)