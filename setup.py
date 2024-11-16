from setuptools import setup, find_packages

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

VERSION = '1.0.7.2' 
DESCRIPTION = 'TnSeeker'
LONG_DESCRIPTION = 'Versatile Python3 module for processing and analyzing anything related with Tn-Seq. Requires bowtie2 to be callable from terminal. File wise, only fastq and genbank annotation files are required.'

setup(
       # the name must match the folder name 
        name="tnseeker", 
        version=VERSION,
        author="Afonso M Bravo",
        author_email="<afonsombravo@hotmail.com>",
        url='https://github.com/afombravo/tnseeker',
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type='text/markdown',
        packages=find_packages(),
        install_requires=["numpy == 1.26.4",
                           "matplotlib == 3.8.2",
                           "pandas ==  2.2.0", 
                           "numba == 0.59.0",
                           "biopython == 1.83",
                           "datetime == 5.4",
                           "argparse == 1.4.0",
                           "regex == 2023.12.25",
                           "multiprocess == 0.70.16",
                           "pathlib == 1.0.1",
                           "scipy == 1.12.0",
                           "seaborn == 0.13.2",
                           "statsmodels == 0.14.1",
                           "colorama"],
        
        entry_points={
        'console_scripts': [
            'tnseeker=tnseeker.__main__:main',  # Replace `2fast2q` with the command name you want to use
            ],
        },

        keywords=['tn-seq', 'essentiality','rb-seq'],
        classifiers= [
            "Development Status :: 3 - Alpha", #change here for when the product is more finished
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3",
            "Operating System :: OS Independent",
        ],
        include_package_data=True,
        package_data={'': ['data/*','data/test/*']},
)