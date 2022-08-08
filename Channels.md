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

**TIP:** Two asterisks, i.e. `**`, works like `*` but crosses directory boundaries. This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns.

For example:

```
files = Channel.fromPath( 'data/**.fa' )
moreFiles = Channel.fromPath( 'data/**/*.fa' )
pairFiles = Channel.fromPath( 'data/file_{1,2}.fq' )
```

The first line returns a channel emitting the files ending with the suffix `.fa` in the data folder and recursively in all its sub-folders. While the second one only emits the files which have the same suffix in any sub-folder in the data path. Finally the last example emits two files: `data/file_1.fq` and `data/file_2.fq`.

**NOTE:** As in Linux Bash, the `*` wildcard does not catch hidden files (i.e. files whose name starts with a `.` character).

In order to include hidden files, you need to start your pattern with a period character or specify the `hidden: true` option. For example:

```
expl1 = Channel.fromPath( '/path/.*' )
expl2 = Channel.fromPath( '/path/.*.fa' )
expl3 = Channel.fromPath( '/path/*', hidden: true )
```

The first example returns all hidden files in the specified path. The second one returns all hidden files ending with the `.fa` suffix. Finally the last example returns all files (hidden and non-hidden) in that path.

By default a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) pattern only looks for regular file paths that match the specified criteria, i.e. it won’t return directory paths.

You may use the parameter `type` specifying the value `file`, `dir` or `any` in order to define what kind of paths you want. For example:

```
myFileChannel = Channel.fromPath( '/path/*b', type: 'dir' )
myFileChannel = Channel.fromPath( '/path/a*', type: 'any' )
```

The first example will return all *directory* paths ending with the `b` suffix, while the second will return any file and directory starting with a `a` prefix.

|Name|Description|
|----|-----------|
|glob | When `true` interprets characters `*`, `?`, `[]` and `{}` as glob wildcards, otherwise handles them as normal characters (default: `true`)|
|type | Type of paths returned, either `file`, `dir` or `any` (default: `file`)|
|hidden | When true includes hidden files in the resulting paths (default: `false`)|
|maxDepth | Maximum number of directory levels to visit (default: *no limit*) |
|followLinks | When `true` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: `true`) |
|relative | When `true` returned paths are relative to the top-most common directory (default: `false`) |
|checkIfExists | When `true` throws an exception of the specified path do not exist in the file system (default: `false`) |

**NOTE:** Multiple paths or glob patterns can be specified using a list:

```
Channel.fromPath( ['/some/path/*.fq', '/other/path/*.fastq'] )
```

### fromFilePairs

The `fromFilePairs` method creates a channel emitting the file pairs matching a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) pattern provided by the user. The matching files are emitted as tuples in which the first element is the grouping key of the matching pair and the second element is the list of files (sorted in lexicographical order). For example:

```
Channel
    .fromFilePairs('/my/data/SRR*_{1,2}.fastq')
    .view()
```

will produce an output similar to the following:

```
[SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
[SRR493367, [/my/data/SRR493367_1.fastq, /my/data/SRR493367_2.fastq]]
[SRR493368, [/my/data/SRR493368_1.fastq, /my/data/SRR493368_2.fastq]]
[SRR493369, [/my/data/SRR493369_1.fastq, /my/data/SRR493369_2.fastq]]
[SRR493370, [/my/data/SRR493370_1.fastq, /my/data/SRR493370_2.fastq]]
[SRR493371, [/my/data/SRR493371_1.fastq, /my/data/SRR493371_2.fastq]]
```

**NOTE:** The glob pattern must contain at least one `*` wildcard character.

Alternatively it is possible to implement a custom file pair grouping strategy providing a closure which, given the current file as parameter, returns the grouping key. For example:

```
Channel
    .fromFilePairs('/some/data/*', size: -1) { file -> file.extension }
    .view { ext, files -> "Files with the extension $ext are $files" }
```

Table of optional parameters available:

|Name|Description|
|----|-----------|
|type|Type of paths returned, either file, dir or any (default: file)|
|hidden|When `true` includes hidden files in the resulting paths (default: `false`)|
|maxDepth|Maximum number of directory levels to visit (default: no *limit*)|
|followLinks|When `true` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: `true`)|
|size|Defines the number of files each emitted item is expected to hold (default: 2). Set to `-1` for any.|
|flat|When `true` the matching files are produced as sole elements in the emitted tuples (default: `false`).|
|checkIfExists|When `true` throws an exception of the specified path do not exist in the file system (default: `false`)|

**NOTE:** Multiple glob patterns can be specified using a list:

```
Channel.fromFilePairs( ['/some/data/SRR*_{1,2}.fastq', '/other/data/QFF*_{1,2}.fastq'] )
```

### fromSRA

### of

The `of` method allows you to create a channel emitting any sequence of values that are specified as the method argument. For example:

```
ch = Channel.of{1, 2, 3, 4}
ch.view{"value: $it"}
```

The first line in this example creates a variable `ch` which holds a channel object. This channel emits the values specified as a parameter in the `of` method. Thus, the second line prints the following:

```
value: 1
value: 3
value: 5
value: 7
```

Ranges of values are expanded accordingly:

```
ch = Channel
    .of(1..22,'X','Y')
    .view()
```

This prints:

```
1
2
3
.
.
.
22
X
Y
```

**NOTE** This feature requires Nextflow version 19.10.0 or later.

See also: [fromList](https://www.nextflow.io/docs/latest/channel.html#fromlist) factory method.

### value

The `value` factory method is used to create a value channel. An optional not `null` argument can be specified to bind the channel to a specific value. For example:

```
expl1 = Channel.value()
expl2 = Channel.value( 'Hello there' )
expl3 = Channel.value( [1,2,3,4,5] )
```

The first line in the example creates an ‘empty’ variable. The second line creates a channel and binds a string to it. Finally, the last one creates a channel and binds a list object to it that will be emitted as a sole emission.

### watchPath

The `watchPath` method watches a folder for one or more files matching a specified pattern. As soon as there is a file that meets the specified condition, it is emitted over the channel that is returned by the `watchPath` method. The condition on files to watch can be specified by using `*` or `?` wildcard characters i.e. by specifying a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) path matching criteria.

For example:

```
Channel
    .watchPath('/path/*.fa')
    .subscribe{println "Fasta file: $it"}
```

By default it watches only for new files created in the specified folder. Optionally, it is possible to provide a second argument that specifies what event(s) to watch. The supported events are:

|Name|Description|
|----|-----------|
|`create`| A new file is created (default)|
|`modify`| A file is modified|
|`delete`| A file is deleted|

You can specify more than one of these events by using a comma separated string as shown below:

```
Channel
    .watchPath( '/path/*.fa', 'create,modify' )
    .subscribe { println "File created or modified: $it" }
```

**WARNING:** The `watchPath` factory waits endlessly for files that match the specified pattern and event(s), which means that it will cause your pipeline to run forever. Consider using the `until` operator to close the channel when a certain condition is met (e.g. receiving a file named `DONE`).

See also: [fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath) factory method.

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

