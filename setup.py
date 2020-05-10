import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="bitclust",
    version="0.0.11",
    description="Fast and memory-efficient clustering of long Molecular Dynamics",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/rglez/bitclust",
    author="Roy González-Alemán",
    author_email="roy.gonzalez-aleman@u-psud.fr",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
    packages=["bitclust"],
    include_package_data=True,
    install_requires=["bitarray>=1.2.1", "numpy", "pandas", "matplotlib"],
    entry_points={
        "console_scripts": [
            "bitclust = bitclust.__main__:main",
        ]
    },
)
