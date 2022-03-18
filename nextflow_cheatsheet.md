# Nextflow DSL2 getting-started cheatsheet

First of all, start using the new syntax `DSL2`. To achieve this, you have to add the following line:

```
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
```

## Key components of a pipeline
- **Inputs**: in the context of Nextflow, they are stored in [Channels](https://www.nextflow.io/docs/latest/channel.html).
- **Data processing steps**: in the context of Nextflow DLS2, these are defined as [Processes](https://www.nextflow.io/docs/latest/process.html) and organized into [Workflows](https://www.nextflow.io/docs/latest/dsl2.html#workflow). 
    - The output from one process is stored into a channel and can be piped into the next process. 
    - If an input channel contains multiple elements, Nextflow will automatically run one process for each in parallel.  
    - Each execution of a process will run in its own working directory, with input often created as symbolic link(symlink) to the original file. 
    - Output files are by default generated in the working directory, and [copied to](https://www.nextflow.io/docs/latest/process.html?highlight=publish#publishdir) the specified directory. 
- **Computing environment and resources** are set up to run the pipeline: 
    - Will it run on a local machine, a high-performance computing environment (HPC) or the cloud?
    - How much CPU or memory will be used? 
In the context of Nextflow, these are specified by [Executors](https://www.nextflow.io/docs/latest/executor.html) and often stored in separate [configuration files](https://www.nextflow.io/docs/latest/config.html). 


### A minimal example
```
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Create input channel. Each .txt file is one element.
input_ch = Channel.fromPath( "*.txt" )  

workflow { 
  input_ch | p1
  p1.out.output_ch2 | p2 | p3
}

process p1 {

  // Run locally instead of HPC or Cloud
  executor 'local'  
  
  input:
    file(x) 

  output:
    file("head.txt")
    /* p2.out will include all output files 
    from p2, whereas emit gives this 
    specific channel a name */
    file("tail.txt"), emit: output_ch2  

  """
  head $x > head.txt
  tail $x > tail.txt
  """
}

process p2 {

  executor 'local'
  // Copy files out of the working directory
  publishDir 'output_folder', mode: 'copy'  

  input:
    file(y)

  output: 
    file("*.gz")

  """
  gzip $y
  """
}
```

### Channels

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

The channel contents can be combined or modified using `operators`. 

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

**Until:** The `until` operator returns all the values until a certain
condition is met (the last value that meets the condition is NOT included)

```
Channel
  .from(1..100)
  .until({it == 49})
  .view()
```
#### Transforming operators

**Transforming operators** are operators that get the items emitted by a channel
and they transform them to new values

**map:** This operator applies a chosen function to every item of a channel

```
Channel
  .from(1, 2, 3, 4)
  .map({it * it})
  .subscribe onNext: {println it}, onComplete: {println "Done!"}
```

**flatMap:** This operator is like map but here instead of a list of items,
each item is returned individually

```
Channel
  .from(1, 2, 3, 4)
  .flatMap({it * it})
  .view()
```

[comment]: # (To be added // reduce)

**groupTuple:** This operator groups items emitted by a channel using a mapping
function which associates a value with a key

```
Channel
  .from( [1, 'A'], [1, 'B'], [2, 'A'], [2, 'c'] )
  .groupTuple()
  .view()
```

**collate:** The collate operator transforms a channel
in such a way that the emitted items are grouped in tuples containing n number
of items where n is specified by the user

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

**Collect:** This operator collect all the items emitted from a channel to a
list and returns them as a single list object

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

**flatten:** This operator transforms a channel so that each item is emitted separately even if it originally belongs to a collection or an array

```
Channel
  .from(1, [3, 4], 8, [34, 35, 36])
  .flatten()
  .view()
```
#### Combining operators

The `combining operators` combine the emitted values from multiple channels

**join:** The `join` operator creates a channel that joins together the items emitted by two channels when a matching key exists

```
ch1 = Channel.from(['X', 1], ['Y', 2])
ch2 = Channel.from(['X', 6], ['Y', 3])
ch1.join(ch2).view()
```

**mix:** The `mix` operator combines the items of more than one channels into one

```
c1 = Channel.from( 1,2,3 )
c2 = Channel.from( 'a','b','c' )
c3 = Channel.from( 'y','z' )
c1.mix(c2, c3).view()
```

**collectFile:** The `collectFile` operator collects the channel emissions and saves them into one or more files

```
Channel
    .from('alpha', 'beta', 'gamma')
    .collectFile(name: 'sample.txt', newLine: true)
    .subscribe {
        println "Entries are saved to file: $it"
        println "File content is: ${it.text}"
    }
```

**combine:** The `combine` operator returns the Cartesian product of items emitted by two channels

```
ch1 = Channel.from(1..5)
ch2 = Channel.from('A'..'C')
ch1.combine(ch2).view()
```

**concat:** The `concat` operator concatenates and returns the items from two or more channels but unlike `mix` it
retains the order

```
a = Channel.from('a','b','c')
b = Channel.from(1,2,3)
c = Channel.from('p','q')
c.concat( b, a ).view()
```

### Processes

Here are the major components of a `process`:

```
process < name > {

  [ directives ]        

  input:                
  < process inputs >

  output:               
  < process outputs >

  when:                 
  < condition >

  [script|shell|exec]:  
  """
  < user script to be executed >
  """
}
```

The `name`, the `input` and `output` of the process are specified. Conditionals (`when`) can also be specified, so that the process runs when certain conditions are met. Note that when using `DSL2`,  there is no need for including the words `from` and `into` for creating the input and output channels.The `script` block defines the command to be executed. This block is interpreted by default as `bash` script but other code can be used too if the `Shebang (#!)` declaration is present:

```
process pyStuff {
  script:
  """
  #!/usr/bin/env python
  print("Hello world!")
  """
}
```

If, instead of using `"""`, `'''` are used, then `Bash` variables can be directly called without escaping `$`. For example:

```
process bar {
  script:
  '''
  echo $PATH | tr ':' '\n'
  '''
}
```

Insted of `script`, `shell` can be used in order to mix `Bash` variables and `Nextflow` variables. In this case, `Nextflow` variables should be defined using the `!{..}` syntax:

```
params.data = 'le monde'

process baz {
  shell:
  '''
  X = 'Bonjour'
  echo $X !{params.data}
  '''
}
```

Here is a simple process which prints the corresponding input to the console:

```
process printWord{
  input:
  val x

  output:
  stdout

  script:
  """
  echo $x
  """
}
```
Here is another one which converts the input to uppercase:

```
process upper{
  input:
  val x

  output:
  stdout

  script:
  """
  echo "$x" | tr '[a-z]' '[A-Z]'
  """
}
```
The output of the first process `hello` becomes input for the second one, converting lowercase `hello` to uppercase `HELLO`.

Files are often used as inputs and/or outputs in processes and thus knowing some file attributes can be extremely useful. Some commonly used file attributes are given in the table below:

**Some useful file attributes**

|Attribute|What it does|
|---------|------------|
| getName | gets the name of the file (ignores the path) |
| getBaseName | gets the file name without its extension |
| getSimpleName | gets the file name after removing any extension |
| getExtension | gets the extension of the file |
| exists | check if the file exists |
| isFile | returns `true` if it is a regular file |
| isDirectory | returns `true` if it is a directory |

### Workflows

`Workflows` are sets of processes that take some inputs through a series of steps in order to produce a certain output. In a workflow, the contents of a channel can become input to another process. Thus, multiple processes can be chained. As an example, in the following workflow, a channel with the word `Hello` is created by the process `printWord` and the content of this channel is passed to the process `upper` to print `HELLO` 

```
workflow {
  a = printWord("hello")
  upper(a).view()
}
```

Having defined `inputs`, `channels` and `workflows`, we can now write a complete Nextflow script and run it. Here is a very simple but complete Nextflow script which writes a greeting message to a file called `hello.txt`:

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

To include `log` information, we include the following:

```
log.info """
         """
         .stripIndent()
```

Between the triple brackets (`"""`), we include the parameters and their usage, as well as the outputs.

To run the script we first save the code snippet above to a file with `.nf` extension (e.g., `main.nf`). After having Nextflow installed and assuming that the `main.nf` file is in our working directory, we run the following on the terminal: `nextflow run main.nf`.

To publish the output text file in a directory, we need to use `publishDir` as shown below:

```
publishDir "Hello", copy: true
```

### Modules in DLS2

A main advantage of the `DSL2` syntax extension is the ability to write and use `modules`. Modules can be included and shared across workflows. Thus, code repetition can be avoided and Nextflow pipelines become more succinct. In addition, the nf-core community maintains high-quality modules for commonly used tools, which can be found here: [Nextflow modules](https://github.com/nf-core/modules/tree/master/modules). 

Modules may contain process, function, and workflow definitions. Components defined in a module, can be imported in another Nextflow script using the keyword `include` as shown in the example below:

```
include { foo } from './some/module' 

workflow {
    data = Channel.fromPath('/data/*.txt')
    foo(data)
}    
```

In this example, a process called `foo` which is present in the module `./some/module` is invoked. The process `foo` takes `data` as input.

When multiple components need to be included from the same `module`, the components can be specified in the same inclusion. Their names need to be separated by `;` as shown below:

```
include { foo; bar } from './some/module' 

workflow {
    data = Channel.fromPath('/data/*.txt')
    foo(data)
    bar(data)
}    
```
### Functions in DSL2

DSL2 allows us to write functions, such as the ones shown below:

```
// Write a function
def print_on_console(x) {
  println x
}

print_on_console("Hello!")
```
The function returns the last evaluated expression, unless a `return` statement is provided explicitly, as in the example given below:

```
def fib (x) {
  if (x <= 1)
    return x
  else
    fib(x - 1) + fib(x - 2)
}

println fib(3)
```

## Sources

- [Official Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)
- [Seqera Labs Nextflow Training](https://training.seqera.io/)
- [Nextflow cheatsheet by Dan Lu](https://github.com/danrlu/Nextflow_cheatsheet)

## Contributors

- Charalampos (Harris) Lazaris 
[GitHub](https://github.com/chlazaris) 
[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/chlazaris.svg?style=social&label=Follow%20%40chlazaris)](https://twitter.com/chlazaris)
- Saba Nafees 
[GitHub](https://github.com/snafees)
[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/saba_nafees314.svg?style=social&label=Follow%20%40saba_nafees314)](https://twitter.com/saba_nafees314)
- Dan Lu
[GitHub](https://github.com/danrlu)

