#!/bin/env nextflow

// Enable DSL2 syntax 
nextflow.enable.dsl=2

samples = Channel.from('10', '100', '1000')
names = Channel.from('name1', 'name2', 'name3')
samples.view()
names.view()

// This is how you combine the two channels
// (equivalent to deprecated spread)
names.
    combine(samples)
    .view()