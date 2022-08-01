# Channels

Nextflow is based on the Dataflow programming model in which processes communicate through channels.

A channel has two major properties:

1. Sending a message is an asynchronous operation which completes immediately, without having to wait for the receiving process.

2. Receiving data is a blocking operation which stops the receiving process until the message has arrived.

## Channel types

Nextflow distinguish two different kinds of channels: *queue channels* and *value channels*.

### Queue channel

A queue channel is a non-blocking unidirectional FIFO queue which connects two processes or operators.

A queue channel is usually created using a factory method such as a [from](https://www.nextflow.io/docs/latest/channel.html#from), [fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath), etc. or chaining it with a channel operator such as [map](https://www.nextflow.io/docs/latest/operator.html#operator-map), [flatMap](https://www.nextflow.io/docs/latest/operator.html#operator-flatmap), etc.

Queue channels are also created by process output declarations using the `into` clause (a necessity in DSL1).

**NOTE:** The definition implies that a queue channel can only be used once as a process output and once as a process input.

If you need to connect a process output channel to more than one process or operator use the into operator to create two (or more) copies of the same channel and use each of them to connect a separate process.

### Value channel

A *value channel* a.k.a. *singleton channel* by definition is bound to a single value and it can be read unlimited times without consuming its content.

**NOTE:** For this reason a value channel can be used as input by more than one process.

A value channel is created using the [value](https://www.nextflow.io/docs/latest/channel.html#value) factory method or by operators returning a single value, such us [first](https://www.nextflow.io/docs/latest/operator.html#operator-first), [last](https://www.nextflow.io/docs/latest/operator.html#operator-last), [collect](https://www.nextflow.io/docs/latest/operator.html#operator-collect), [count](https://www.nextflow.io/docs/latest/operator.html#operator-count), [min](https://www.nextflow.io/docs/latest/operator.html#operator-min), [max](https://www.nextflow.io/docs/latest/operator.html#operator-max), [reduce](https://www.nextflow.io/docs/latest/operator.html#operator-reduce), [sum](https://www.nextflow.io/docs/latest/operator.html#operator-sum).

**NOTE:** A value channel is implicitly created by a process when an input specifies a simple value in the `from` clause. Moreover, a value channel is also implicitly created as output for a process whose inputs are only value channels.

For example

```
// Enable DSL2
nextflow.enable.dsl=2

process foo {
    input:
    val x

    output:
    file 'x.txt'

    shell:
    """
    echo !{x} > x.txt
    """
}

workflow {
    foo('hello')
}
```

The process in the above snippet declares a single input which implicitly is a value channel.

See also: [Understand how multiple input channels work](https://www.nextflow.io/docs/latest/process.html#process-understand-how-multiple-input-channels-work).

## Channel factory

### create

### empty

### from

### fromPath

### fromFilePairs

### fromSRA

### of

### value

### watchPath

## Binding values

## Observing events

