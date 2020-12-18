import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Ekman",
    version="1.0",
    author="Ueslei Sutil",
    author_email="ueslei@outlook.com",
    description="Python toolbox to handle with models outputs.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/uesleisutil/Ekman",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.5',
    install_requires = [
        "numpy>=1.19.2",
        "dask==0.15.2",
        "netCDF4==1.5.3",
        "wrf-python==1.3.2",
        "matplotlib == 3.3.2"
    ]
)