import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="cellmarkerpipe",
    version='0.0.0',
    description="A pipeline to select marker genes.",
    author="Yinglu Jia",
    packages=find_packages(include=['pipeline']),
    long_description=read('README.md'),
    entry_points = {'console_scripts': ['cellMarkerPipe=pipeline.cellMarkerPipe:main']},
    install_requires=[
        'pandas',
        'numpy'
        ]
)

