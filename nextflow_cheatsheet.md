# Notes for Nextflow - nf-core Hackathon March 2022

## Here are some notes on how to use Nextflow

First of all, start using the new syntax `DSL2`. To achieve this, you have to
use the following line:

```
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

```

Now you can write a simple process and the corresponding output:

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
  cat $x > hello.txt
}

// Specify the workflow
workflow {
    writeText(greeting)
}
```
