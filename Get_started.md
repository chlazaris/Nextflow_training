# Get started

## Requirements

Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and [Java 11 (or later, up to 18)](https://www.oracle.com/java/technologies/downloads/) to be installed.

For the execution in a cluster of computers, the use of a shared file system is required to allow the sharing of tasks input/output files.

Nextflow can also be run on Windows through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux)

## Installation

Nextflow is distributed as a self-installing package, which means that it does not require any special installation procedure.

It only needs two easy steps:

1. Download the executable package by copying and pasting either one of the following commands in your terminal window: `wget -qO- https://get.nextflow.io | bash`.

Or, if you prefer *curl*: `curl -s https://get.nextflow.io | bash`

This will create the `nextflow` main executable file in the current directory.

2. Make the binary executable on your system by running: `chmod +x nextflow`

An alternative would be to install Nextflow using [`conda`](https://anaconda.org/bioconda/nextflow) or [`mamba`](https://gist.github.com/ckandoth/982ce140b4dd9d6bf72a780c05a549a3)

Optionally, move the `nextflow` file to a directory accessible by your `$PATH` variable (this is only required to avoid remembering and typing the full path to `nextflow` each time you need to run it).

## Updates

Having Nextflow installed in your computer you can update to the latest version using the following command:

`nextflow self-update`

## Stable & Edge releases

A *stable* version of Nextflow is released on a six-months basic schedule, in the 1st and 3rd quarter of every year.

Along with the stable release, an `edge` version is released on a monthly basis. This version is useful to test and use most recent updates and experimental features.

To use the latest `edge` release, run the following snippet in your shell terminal:

```
export NXF_EDGE=1
nextflow self-update
```

## Your first script

### Modify and resume

### Pipeline parameters
