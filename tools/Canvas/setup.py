# Setup Python module
from setuptools import setup, find_packages

modules = ["Canvas." + p for p in sorted(find_packages("./Canvas"))]

setup(
    name="Canvas",
    version="0.1",
    description="Genome Visualisation Tool",
    url="https://github.com/breons/canvas",
    author="Breon Schmidt",
    license="MIT",
    packages=["Canvas", *modules],
    zip_safe=False,
    include_package_data=True,
    install_requires=[
		"kaleido==0.1.0",
        "numpy==1.18.1",
        "pandas==1.2.3",
		"plotly==4.14.3",
		"pybedtools==0.8.1",
		"pysam==0.15.4"
		# from functools import reduce
    ]
)