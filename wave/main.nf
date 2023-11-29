#!/usr/bin/env nextflow

/* This is a Nextflow pipeline that
uses cowpy to print a message to the screen */

// Activate DSL2 
nextflow.enable.dsl=2

// Specify the cowpy process
process Cowpy {
    // Specify the container
    conda "/Users/lazarch2/Desktop/wave/cowpy_env.yaml"

    output:
    stdout

    script:
    """
    python cowpy.py
    """
}

workflow {
    Cowpy()
}
