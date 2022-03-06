# Notes for Nextflow - nf-core Hackathon March 2022

## Here are some notes on how to use Nextflow

First of all, start using the new syntax `DSL2`. To achieve this, you have to
use the following line:

```
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

```

Now you can write a simple process and the corresponding output: s

```
params.greetings="Hello world"
greeting = Channel.from(params.greetings)

// Write the processes
process writeText {
  input:
  val x

  output:
  file "hello.txt"

  script:
  """
  echo ${x} > hello.txt
  """
}

// Specify the workflow
workflow {
    writeText(greeting)
}
```

Now, order to run the workflow after having Nextflow installed, you need
to run the following: `nextflow run main.nf` where `main.nf` is the name
of the script containing the code above.

To publish the output text file in a directory, we need to use `publishDir`
as shown below:

```
publishDir "Hello", copy: true
```

**Some useful file attributes**

|Attribute|What it does|
|---------|------------|
| getName | gets the name of the file (ignores the path) |
| getBaseName | gets the file name without its extension |
| getSimpleName | gets the file name after removing any extension |
| getExtension | gets the extension of the file |
| exists | check if the file exists |
| isFile | returns `true` if it a regular file |
| isDirectory | returns `true` if it is directory |

To include `log` information for the pipeline, we include the following in the
`.nf` file:

```
log.info """
         """
         .stripIndent()
```

Between the triple brackets (`"""`), we include the parameters and their usage,
as well as the outputs.
