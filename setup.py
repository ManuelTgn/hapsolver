from setuptools import setup, find_packages

# Reading the long description from the README file (optional)
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="hapsolver",  
    version="0.1.0",  # Initial release version
    author="ManuelTgn",  
    author_email="manu.tognon@gmail.com",  
    description="Haplotype reconstruction tool",  # A brief description of the package
    long_description=long_description,  # Detailed description from README
    long_description_content_type="text/markdown",  # Format of the long description (Markdown here)
    url="https://github.com/ManuelTgn/hapsolver",  
    packages=find_packages(),  # Automatically find all packages in the current directory
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
    ],
    python_requires='>=3.8',  
    install_requires=[  
        "pysam",
    ],
    extras_require={ 
        "dev": [
            "pytest",
            "black",
        ],
        "docs": [
            "sphinx",
            "sphinx_rtd_theme",
        ],
    },
    include_package_data=True,  
    zip_safe=False,  
)
