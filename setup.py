from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pydemult',
      version='0.6',
      description='Streamed and parallel demultiplexing of fastq files in python',
      long_description=readme(),
      url='https://github.com/jenzopr/pydemult',
      author='Jens Preussner',
      author_email='jens.preussner@mpi-bn.mpg.de',
      license='MIT',
      packages=['pydemult'],
      entry_points = {
        'console_scripts': ['pydemult = pydemult.pydemult:demultiplex',
                            'HTOcount = pydemult.HTO:count']
      },
      scripts = [],
      setup_requires=["numpy","Cython"],
      install_requires=[
        'wheel',
        'regex',
        'pandas',
        'mputil',
        'scvi-tools',
        'anndata',
        'python-Levenshtein'
      ],
      classifiers = [
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3'
      ],
      zip_safe=False,
      include_package_data=True)
