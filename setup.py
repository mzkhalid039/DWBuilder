from setuptools import setup, find_packages

setup(
    name='dwbuilder',
    version='2.0.0',
    author='M.Z.Khalid, S.M.Selbach',
    author_email='zeeshan.khalid039@gmail.com',
    description='A package for developing domain wall structures for atomistic calculations',
    license='MIT',
    url='https://github.com/mzkhalid039/DWBuilder',
    packages=find_packages(include=['scripts']),
    py_modules=["main"],  # Ensure main.py is included
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        'ase',
        'numpy<2',
        'matplotlib',
        'pymatgen',
        'colorama'
    ],
    entry_points={
        'console_scripts': [
            'dwbuilder = main:main',
        ]
    }
)
