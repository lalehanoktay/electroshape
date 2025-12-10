## Introduction

ElectroShape is a ligand-based virtual screening method that combines shape, electrostatic information, and lipophilicity into a single, unified framework. It builds upon the ultra-fast shape recognition (USR) method by extending it to higher dimensions.

The original **ElectroShape** method, introduced by **Armstrong et al. in 2010**, represents molecules in a 4-dimensional space, using the three spatial coordinates (x, y, z) and partial charge (q) as the fourth dimension. This allows for a more descriptive representation of molecules, capturing both their shape and electrostatic properties.

In **2011**, the method was further improved by incorporating atomic lipophilicity (alogP) as a **fifth dimension. This 5D approach was shown to improve the accuracy of virtual screening. This tool implements the 5D ElectroShape method.

## Citations

1.  Armstrong, M. S., Morris, G. M., Finn, P. W., Sharma, R., Moretti, L., Cooper, R. I., & Richards, W. G. (2010). ElectroShape: fast molecular similarity calculations incorporating shape, chirality and electrostatics. *Journal of Computer-Aided Molecular Design, 24*(9), 789–801.
2.  Armstrong, M. S., Finn, P. W., Morris, G. M., & Richards, W. G. (2011). Improving the accuracy of ultrafast ligand-based screening: incorporating lipophilicity into ElectroShape as an extra dimension. *Journal of Computer-Aided Molecular Design, 25*(8), 785–790.

## Installation

To install the package, clone the repository and install the dependencies:

```bash
git clone <repository-url>
cd electroshape
pip install -r requirements.txt
pip install .
```

## Usage

The tool can be run from the command line:

```bash
electroshape --sdf <input.sdf> --out <output.csv>
```

### Arguments

-   `--sdf`: Path to the input SDF file (can be gzipped, `.sdf.gz`).
-   `--out`: Path to the output CSV file (can be gzipped, `.csv.gz`).
-   `--workers`: Number of CPU cores to use (default: number of cores - 1).
-   `--chunk`: Number of molecules to process in a batch (default: 5000).
-   `--charge`: Method for partial charge calculation (`mmff` or `gasteiger`, default: `mmff`).
