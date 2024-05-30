from setuptools import setup, find_packages
from os.path import abspath, dirname

setup(
    name='dwbuilder',
    version='1.1.0',
    authors = 'M.Z.Khalid, S.M.Selbach',
    author_email='zeeshan.khalid039@gmail.com',
    description='A package for developing domain wall structures for atomistic calculations',
    license='MIT',
    url='https://github.com/mzkhalid039/DWBuilder',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=['ase','numpy','matplotlib', 'pymatgen', 'colorama'],
    entry_points={
        'console_scripts': [
            'dwbuilder = scripts.dwbuilder:main',
            'dbuilder = scripts.dbuilder:main',
            'hibuilder = scripts.hibuilder:main',
            'polarization = scripts.polarization:main',
        ]
    }
)

