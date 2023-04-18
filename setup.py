from setuptools import setup, find_packages
from os.path import abspath, dirname

setup(
    name='DWBuilder',
    version='1.0.0',
    author='M.Z.Khalid',
    author_email='muhammad.z.khalid@ntnu.no',
    description='A package for developing domain wall structures for atomistic calculations',
    license='MIT',
    packages=find_packages(),
    install_requires=['ase','numpy','matplotlib'],
    entry_points={
        'console_scripts': [
            'dwbuilder = scripts.dwbuilder:main',
            'dbuilder = scripts.dbuilder:main',
        ]
    }
)

