
# MultiQC for Chemical GlycoBiology
### Aggregate chemical glycobiology results across many samples into a single report.


MultiQC is a tool to create a single report with interactive plots
for multiple bioinformatics analyses across many samples.

This fork includes [computational chemistry, glycoinformatics and chemical biology tools, more documentation here](README-chemicalglycobiology.md)


This fork of MultiQC is written in Python (tested with v 3.6). It is not yet available through pypi or conda.

Currently, supported computational, glycoinformatics and chemical biology tools include:

|Molecular Dynamics               | Quantum Chemistry      | Post-dynamics analysis (e.g. Ring Pucker) | Other |
|---------------------------------|------------------------|---------------------------|----------------------|
|                                 |comp_qm                 |comp_tesselate             |                      |


## Installation

### Prerequisite dev environment
Install Ananconda python3, then create an environment for MultiQC

conda create -n multiqc
source activate multiqc

### clone the code
git clone
pushd

### compile the code
python setup.py install > log.log

## Usage
Once installed, you can use MultiQC by navigating to your analysis directory
(or a parent directory) and running the tool:
```bash
multiqc .
```

Select log files from the data dir, only use the quantum analysis module and force creation of a report (overwrites existing)
```bash
multiqc data/*.log -m comp_qm -f
```


## Development

### New modules
Create a new module, e.g. add module newmodule in multiqcc say multiqc/modules/newmodule/*
A great example to work from is kallisto

Then include newmodule in:
 - multiqc/utils/search_patterns.yaml
 - setup.py


# Future ideas (which may be fixed by core multiqc or highcharts dev)

- render windows separately for large report. i.e. turn pieces on and off for memory intensive report
- fix colour scales to include negative numbers
