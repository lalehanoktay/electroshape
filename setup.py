from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="electroshape",
    version="1.0.0",
    author="<Your Name>",
    author_email="<your.email@example.com>",
    description="A tool to generate 5D ElectroShape descriptors for molecules.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="<your-github-repo-url>",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
        "numpy<2",
        "rdkit-pypi",
        "tqdm",
    ],
    entry_points={
        "console_scripts": [
            "electroshape=electroshape.cli:main",
        ],
    },
)
