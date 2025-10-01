from setuptools import setup, find_packages

setup(
    name="mldmtd",
    version="1.1.3",
    author="Emmanuel Romero",
    author_email="romeroqe@gmail.com",
    description="Methodology to locate the minimum and maximum depth of the thermocline",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/romeroqe/mld-mtd",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "gsw",
        "netCDF4",
        "scipy",
        "scikit-learn"
    ],
    license="CC-BY-4.0",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
    ],
)