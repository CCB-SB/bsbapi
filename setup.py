from setuptools import setup

setup(
    name='bsbapi',
    version='0.1.0',
    author='Georges Schmartz',
    packages=['bsbapi'],
    install_requires=['requests', 'pandas'],
    licence='LICENSE.txt',
    description="Package for API usage of the Busybee webserver"
    )
