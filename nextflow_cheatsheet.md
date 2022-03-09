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

### Specifying the inputs

A way to specify a channel from different values, is the following:

```
value_ch = Channel.from(1,2,3)
```

Another way to achieve the same result is:

```
Channel.from(1,2,3)
  .set{value_ch}
```

To create a `value` channel, we use the `value` factory method. For example:

```
example_1 = Channel.value()
example_2 = Channel.value('Hello there!')
example_3 = Channel.value([1,2,3,4,5])
```

To create a channel which emits list elements, we can use the `fromList` method:

```
example_1 = Channel.fromList([1,2,3,4])
```

To create a channel from paths, we can use the `fromPath` method:

```
example_1 = Channel.fromPath('/path/file.txt')
```

To check if the file exists, we need to add `checkIfExists: true` as shown
below:

```
example_1 = Channel.fromPath('/path/file.txt', checkIfExists: true)
```

To get the file pairs matching a glob pattern, we need to use the `fromFilePairs`
method:

```
example_1 = Channel.fromFilePairs('/path/*_{1,2}.fastq')
```

Finally, to retrieve records directly from SRA, we use the method `fromSRA`:

```
example_1 = Channel.fromSRA('SRP043510')
```
### Operators

#### Filtering operators

**Filtering operators** are operators that allow to get the emitted elements
from a channel which satisfy certain conditions

**Filter:** The `filter` operator allows to filter results based on a certain
pattern or condition

```
Channel
  .from('a', 'b', 'c', 'aa', 'ab')
  .filter( ~/^a.*/ )
  .view()
```

```
Channel
  .from('a', 'b', 'c', '1', 1, 2, 2.35)
  .filter (Number)
  .view()
```

**Unique:** The `unique` operator allows to return the unique
values from a channel

```
Channel
  .from(1, 2, 3, 1, 4, 'a', 'b', 'a')
  .unique()
  .view()
```

**Distinct:** The `distinct` operator allows to return unique
consecutive values from a channel

```
Channel
  .from(1, 2, 3, 1, 1, 4, 'a', 'b', 'a')
  .distinct()
  .view()
```

**Take:** The `take` operator returns the first n items
emitted by a channel

```
Channel
  .from(1..100)
  .take(10)
  .view()
```

**First:** The `first` operator either returns the first item
or the first one that meets a certain condition

```
Channel
  .from(1, 2, 5, 8)
  .first({it > 4})
  .view()
```

**Last:** The `last` operator returns the last item
of a channel

```
Channel
  .from(1, 2, 5, 8)
  .last()
  .view()
```

**Until:** The `until` operator returns all the values till a certain
condition is met (the last value that meets the condition is NOT included)

```
Channel
  .from(1..100)
  .until({it == 49})
  .view()
```
#### Transforming operators

**Transforming operators** are operators that get the items emitted by a channel and they transform them to new values

**map:** This operator applies a chosen function to every item of a channel

```
Channel
  .from(1, 2, 3, 4)
  .map({it * it})
  .subscribe onNext: {println it}, onComplete: {println "Done!"}
```

**flatMap:** This operator is like map but here instead of a list of items, each item is returned individually

```
Channel
  .from(1, 2, 3, 4)
  .flatMap({it * it})
  .view()
```

[comment]: # (To be added // reduce)

**groupTuple:** This operator groups items emitted by a channel using a mapping function which associates a value
with a key

```
Channel
  .from( [1, 'A'], [1, 'B'], [2, 'A'], [2, 'c'] )
  .groupTuple()
  .view()
```

**collate:** The collate operator transforms a channel
in such a way that the emitted items are grouped in tuples containing n number of items where n is specified by the user

```
Channel
  .from(1..7)
  .collate(3)
  .view()
```
Now, if we want to get rid of the remaining item

```
Channel
  .from(1..7)
  .collate(3, false)
  .view()
```

**buffer:** This is an operator that buffers (subsets)
the values to be returned based on certain conditions

```
// Specify end condition
Channel
  .from(1..100)
  .buffer(5)
  .view()
```

```
// Specify start and end condition
Channel
  .from(1..100)
  .buffer(10, 20)
  .view()
```

```
// Specify size
Channel
  .from(1..100)
  .buffer(size: 3, remainder: false)
  .view()
```

**Collect:** This operator collect all the items emitted from a channel to a list and returns them as a single list object

```
Channel
  .from(1..10)
  .collect()
  .view()
```

**toList:** This operator does what collect does

```
Channel
  .from(1..10)
  .toList()
  .view()
```

**toSortedList:** This operator returns the items in a sorted list

```
Channel
  .from(1, 2, 8, 5, 3, 4)
  .toSortedList()
  .view()
```

**flatten:** This operator transforms a channel so that each item is emitted seperately even if it originally belongs to a collection or an array

```
Channel
  .from(1, [3, 4], 8, [34, 35, 36])
  .flatten()
  .view()
```
