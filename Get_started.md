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

Copy the following example into your favourite text editor and save it to a file named `tutorial.nf`

```
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.str = 'Hello world!'

process splitLetters {
  output:
    path 'chunk_*'

  """
  printf '${params.str}' | split -b 6 - chunk_
  """
}

process convertToUpper {
  input:
    file x
  output:
    stdout

  """
  cat $x | tr '[a-z]' '[A-Z]'
  """
}

workflow {
  splitLetters | flatten | convertToUpper | view { it.trim() }
}
```

This script defines two processes. The first splits a string into 6-character chunks, writing each one to a file with the prefix `chunk_`, and the second receives these files and transforms their contents to uppercase letters. The resulting strings are emitted on the result channel and the final output is printed by the `view` operator.

Execute the script by entering the following command in your terminal:

```
nextflow run tutorial.nf
```

It will output something similar to the text shown below:

```
N E X T F L O W  ~  version 19.04.0
executor >  local (3)
[69/c8ea4a] process > splitLetters   [100%] 1 of 1 ✔
[84/c8b7f1] process > convertToUpper [100%] 2 of 2 ✔
HELLO
WORLD!
```
You can see that the first process is executed once, and the second twice. Finally the result string is printed.

It’s worth noting that the process `convertToUpper` is executed in parallel, so there’s no guarantee that the instance processing the first split (the chunk Hello) will be executed before the one processing the second split (the chunk world!).

Thus, it is perfectly possible that you will get the final result printed out in a different order:

```
HELLO
WORLD
```

### Modify and resume

*Nextflow* keeps track of all the processes executed in your pipeline. If you modify some parts of your script, only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

For the sake of this tutorial, modify the `convertToUpper` process in the previous example, replacing the process script with the string `rev $x`, so that the process looks like this:

```
process convertToUpper {
    input:
    file x from letters

    output:
    stdout result

    """
    rev $x
    """
}
```

Then save the file with the same name, and execute it by adding the `-resume` option to the command line:

```
nextflow run tutorial.nf -resume
```

It will print output similar to this:

```
N E X T F L O W  ~  version 19.04.0
executor >  local (2)
[69/c8ea4a] process > splitLetters   [100%] 1 of 1, cached: 1 ✔
[d0/e94f07] process > convertToUpper [100%] 2 of 2 ✔
olleH
!dlrow
```

You will see that the execution of the process `splitLetters` is actually skipped (the process ID is the same), and its results are retrieved from the cache. The second process is executed as expected, printing the reversed strings.

### Pipeline parameters

Pipeline parameters are simply declared by prepending to a variable name the prefix `params`, separated by dot character. Their value can be specified on the command line by prefixing the parameter name with a double dash character, i.e. `--paramName`

For the sake of this tutorial, you can try to execute the previous example specifying a different input string parameter, as shown below:

```
nextflow run tutorial.nf --str 'Bonjour le monde'
```
The string specified on the command line will override the default value of the parameter. The output will look like this:

```
N E X T F L O W  ~  version 19.04.0
executor >  local (4)
[8b/16e7d7] process > splitLetters   [100%] 1 of 1 ✔
[eb/729772] process > convertToUpper [100%] 3 of 3 ✔
m el r
edno
uojnoB
```
