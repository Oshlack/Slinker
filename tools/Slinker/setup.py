# Setup Python module
from setuptools import setup, find_packages
modules = ["Slinker." + p for p in sorted(find_packages("./Slinker"))]

setup(
    name="Slinker",
    version="0.1b",
    description="Novel Splice Junction Visualisation",
    url="https://github.com/oshlack/slinker",
    author="Breon Schmidt",
    license="MIT",
    packages=["Slinker", *modules],
    zip_safe=False,
    include_package_data=True,
    install_requires=[
        "numpy==1.18.1",
        "pandas==1.2.3",
		"pybedtools==0.8.1",
		"pysam==0.15.4", 
        'Canvas'
    ]
)