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
    install_requires=[
        'pandas',
        'numpy',
        'subprocess'
        ]
)

