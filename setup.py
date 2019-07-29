#!/usr/bin/env python

from setuptools import setup, find_packages


def main():
    with open('README.md') as f:
        readme = f.read()

    with open('LICENSE') as f:
        license = f.read()

    metadata = dict(
        name='Open-Burnup',
        version='0.1.0',
        description='Open Source depletion code',
        long_description=readme,
        author='Julien de Troullioud de Lanversin',
        author_email='jdtdl@princeton.edu',
        url='http://jdtdl.mycpanel.princeton.edu/',
        license=license,
        packages=find_packages(exclude=('docs', 'graphs', 'notes', 'papers', 'test', 'test_nouveau')),
    )

    setup(**metadata)


if __name__ == '__main__':
    main()

