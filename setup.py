from setuptools import setup
import os
import re


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'docko/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def read_author():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'docko/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__author__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def read_email():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'docko/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__author_email__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def read_git():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'docko/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__url__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='docko',
      version=read_version(),
      description='',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author=read_author(),
      author_email=read_email(),
      url=read_git(),
      license='GPL3',
      project_urls={
          "Bug Tracker": read_git(),
          "Documentation": read_git(),
          "Source Code": read_git()},
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords=['docking', 'protein-engineering'],
      packages=['docko'],
      include_package_data=True,
      package_data={
          'docko.deps': [
              'prepare_gpf.py',
              'mgltools_x86_64Linux2_1.5.7',
              'prepare_dpf4.py',
          ],
      },
      entry_points={
          'console_scripts': [
              'docko = docko.__main__:main'
          ]
      },
      install_requires=['babel==2.16.0',
                        'biopython',
                        'h5py',
                        'jupyterlab',
                        'matplotlib',
                        'networkx',
                        'numpy',
                        'opencv-python',
                        'OpenMM==8.1.1',
                        'pandas',
                        'primer3-py==2.0.3',
                        'prompt_toolkit==3.0.47',
                        'protobuf==5.28.1',
                        'pysam==0.22.1',
                        'pytz==2024.1',
                        'PyYAML==6.0.2',
                        'pyzmq==26.2.0',
                        'rdkit',
                        'regex==2024.9.11',
                        'scikit-learn',
                        'sciutil==1.0.3',
                        'sciviso==1.0.9',
                        'seaborn==0.13.2',
                        'statannot==0.2.3',
                        'statsmodels==0.14.2',
                        'tqdm==4.66.5',
                        'sciutil', 
                        'boltz>=2.0.0'],
      python_requires='>=3.8',
      data_files=[("", ["LICENSE"])]
      )
