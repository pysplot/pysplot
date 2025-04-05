from setuptools import setup, find_packages

setup(
    name='pysplot',                      # Package name
    version='0.1.0',                     # Package version
    author='Brent Smith',                 # Author name
    author_email='brent.smith@jhuapl.edu', # Author email
    description='A spatial and temporal plotting package for heliophysics data',
    long_description=open('README.md').read(),  # Long description from README file
    long_description_content_type='text/markdown',
    url='https://github.com/pysplot/pysplot',  # URL to your package's repository
    packages=find_packages(),            # Automatically find packages in your project
    classifiers=[                        # Classifiers for PyPI (useful for search)
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    install_requires=[                   # List of dependencies for the main package
        'sunpy>=5.1.5',                   # SunPy (for solar physics coordinates and data)
        "pandas>=2.2.3",                  # Pandas (for data manipulation)
        "matplotlib>=3.9.0"               # matplotlib for plotting
    ],
    extras_require={                     # Extra dependencies for specific environments
        'whpi': [                         
            "imageio>=2.9.0",
        ],
        'examples': [
            "cdasws>=1.7.42",
            "cdflib>=0.4.9",
            "xarray>=0.20.1",
        ],
        'test': [                        
            'pytest',                    
        ],
    },
    python_requires='>=3.9',              # Python version requirements
)