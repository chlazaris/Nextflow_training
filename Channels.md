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

Channels may be created implicitly by the process output(s) declaration or explicitly using the following channel factory methods.

The available factory methods are:

* [create](https://www.nextflow.io/docs/latest/channel.html#create) (No longer available in DSL2 syntax)
* [empty](https://www.nextflow.io/docs/latest/channel.html#empty)
* [from](https://www.nextflow.io/docs/latest/channel.html#from) (Deprecated-only to be used for compatibility with legacy code)
* [fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath)
* [fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)
* [fromSRA](https://www.nextflow.io/docs/latest/channel.html#fromsra)
* [of](https://www.nextflow.io/docs/latest/channel.html#of)
* [value](https://www.nextflow.io/docs/latest/channel.html#value)
* [watchPath](https://www.nextflow.io/docs/latest/channel.html#watchpath)
* fromList
  
### create

Creates a new channel by using the `create` method, as shown below:

```
channelObj = Channel.create()
```

**WARNING**: Not available in DSL2

### empty

The empty factory method, by definition, creates a channel that doesn’t emit any value.

See also: [ifEmpty](https://www.nextflow.io/docs/latest/operator.html#operator-ifempty) and [close](https://www.nextflow.io/docs/latest/operator.html#operator-close) operators.

### from

**WARNING:** This method is deprecated and should only be used for backward compatibility in legacy code. Use [of](https://www.nextflow.io/docs/latest/channel.html#of) or [fromList](https://www.nextflow.io/docs/latest/channel.html#fromlist) instead.

The `from` method allows you to create a channel emitting any sequence of values that are specified as the method argument, for example:

```
ch = Channel.from( 1, 3, 5, 7 )
ch.subscribe { println "value: $it" }
```

The first line in this example creates a variable `ch` which holds a channel object. This channel emits the values specified as a parameter in the `from` method. Thus the second line will print the following:

```
value: 1
value: 3
value: 5
value: 7
```

The following example shows how to create a channel from a range of numbers or strings:

```
zeroToNine = Channel.from( 0..9 )
strings = Channel.from( 'A'..'Z' )
```

**NOTE:** When the `from` argument is an object implementing the (Java) [Collection](http://docs.oracle.com/javase/7/docs/api/java/util/Collection.html) interface, the resulting channel emits the collection entries as individual items.

Thus, the following two declarations produce an identical result even tough in the first case the items are specified as multiple arguments while in the second case as a single list object argument:

```
Channel.from( 1, 3, 5, 7, 9 )
Channel.from( [1, 3, 5, 7, 9] )
```

But when more than one argument is provided, they are always managed as single emissions. Thus, the following example creates a channel emitting three entries each of which is a list containing two elements:

```
Channel.from( [1, 2], [5,6], [7,9] )
```

### fromPath

You can create a channel emitting one or more file paths by using the fromPath method and specifying a path string as an argument. For example:

```
myFileChannel = Channel.fromPath( '/data/some/bigfile.txt' )
```

The above line creates a channel and binds it to a [Path](http://docs.oracle.com/javase/7/docs/api/java/nio/file/Path.html) object for the specified file.

**NOTE:** `fromPath` does not check whether the file exists (by default).

Whenever the fromPath argument contains a `*` or `?` wildcard character it is interpreted as a glob path matcher. For example:

```
myFileChannel = Channel.fromPath( '/data/big/*.txt' )
```

This example creates a channel and emits as many `Path` items as there are files with `txt` extension in the `/data/big` folder.




### fromFilePairs

### fromSRA

### of

### value

The `value` factory method is used to create a value channel. An optional not `null` argument can be specified to bind the channel to a specific value. For example:

```
expl1 = Channel.value()
expl2 = Channel.value( 'Hello there' )
expl3 = Channel.value( [1,2,3,4,5] )
```

The first line in the example creates an ‘empty’ variable. The second line creates a channel and binds a string to it. Finally, the last one creates a channel and binds a list object to it that will be emitted as a sole emission.

### watchPath

### fromList

The `fromList` method allows you to create a channel emitting the values provided as a list of elements, for example:

```
Channel
    .fromList( ['a', 'b', 'c', 'd'] )
    .view { "value: $it" }
```

prints:

```
value: a
value: b
value: c
value: d
```

See also: [of](https://www.nextflow.io/docs/latest/channel.html#of) factory method.

**NOTE:** This feature requires Nextflow version 19.10.0 or later.

## Binding values

## Observing events

