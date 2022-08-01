# Processes

In Nextflow a *process* is the basic processing primitive to execute a user script.

The process definition starts with the keyword `process`, followed by process name and finally the process *body* delimited by curly brackets. The process body must contain a string which represents the command or, more generally, a script that is executed by it. A basic process looks like the following example:

```
process sayHello {
  """
  echo 'Hello world!' > file
  """
}
```

A process may contain five definition blocks, respectively: directives, inputs, outputs, when clause and finally the process script. The syntax is defined as follows:

```
process < name > {

  [directives]

  input:
    < process input(s) >

  output:
    < process output(s) >

  when:
    < condition >

  [script|shell|exec]:
    < user script to be executed >
}
```

## Script

The `script` block is a string statement that defines the command that is executed by the process to carry out its task.

A process contains one and only one script block, and it must be the last statement when the process contains input and output declarations.

The entered string is executed as a [Bash](http://en.wikipedia.org/wiki/Bash_(Unix_shell)) script in the host system. It can be any command, script or combination of them, that you would normally use in terminal shell or in a common Bash script.

The only limitation to the commands that can be used in the script statement is given by the availability of those programs in the target execution system.

The script block can be a simple string or multi-line string. The latter simplifies the writing of non trivial scripts composed by multiple commands spanning over multiple lines. For example:

```
process doMoreThings {
  """
  blastp -db $db -query query.fa -outfmt 6 > blast_result
  cat blast_result | head -n 10 | cut -f 2 > top_hits
  blastdbcmd -db $db -entry_batch top_hits > sequences
  """
}
```
As explained in the script tutorial section, strings can be defined by using a single-quote or a double-quote, and multi-line strings are defined by three single-quote or three double-quote characters.

There is a subtle but important difference between them. Like in Bash, strings delimited by a `"` character support variable substitutions, whereas strings delimited by `'` do not.

In the above code fragment the `db` variable is replaced by the actual value defined somewhere in the pipeline script.

**WARNING:** Since Nextflow uses the same Bash syntax for variable substitutions in strings, you must manage them carefully depending on whether you want to evaluate a *Nextflow* variable or a *Bash* variable.

When you need to access a system environment variable in your script you have two options. The first choice is as easy as defining your script block by using a single-quote string. For example:

```
process printPath {
  '''
  echo The path is: $PATH
  '''
}
```

The drawback of this solution is that you will not be able to access variables defined in the pipeline script context, in your script block.

To fix this, define your script by using a double-quote string and escape the system environment variables by prefixing them with a back-slash `\` character, as shown in the following example:

```
process doOtherThings {
  """
  blastp -db \$DB -query query.fa -outfmt 6 > blast_result
  cat blast_result | head -n $MAX | cut -f 2 > top_hits
  blastdbcmd -db \$DB -entry_batch top_hits > sequences
  """
}
```

In this example the `$MAX` variable has to be defined somewhere before, in the pipeline script. Nextflow replaces it with the actual value before executing the script. Instead, the `$DB` variable must exist in the script execution environment and the Bash interpreter will replace it with the actual value.

**TIP:** Alternatively you can use the [Shell](https://www.nextflow.io/docs/latest/process.html#process-shell) block definition which allows a script to contain both Bash and Nextflow variables without having to escape the first.

### Scripts a la carte

The process script is interpreted by Nextflow as a Bash script by default, but you are not limited to it.

You can use your favourite scripting language (e.g. Perl, Python, Ruby, R, etc), or even mix them in the same pipeline.

A pipeline may be composed by processes that execute very different tasks. Using Nextflow you can choose the scripting language that better fits the task carried out by a specified process. For example for some processes R could be more useful than Perl, in other you may need to use Python because it provides better access to a library or an API, etc.

To use a scripting other than Bash, simply start your process script with the corresponding [shebang](http://en.wikipedia.org/wiki/Shebang_(Unix)) declaration. For example:

```
process perlStuff {
    """
    #!/usr/bin/perl

    print 'Hi there!' . '\n';
    """
}

process pythonStuff {
    """
    #!/usr/bin/python

    x = 'Hello'
    y = 'world!'
    print "%s - %s" % (x,y)
    """
}
```

### Conditional Scripts

Complex process scripts may need to evaluate conditions on the input parameters or use traditional flow control statements (i.e. `if`, `switch`, etc) in order to execute specific script commands, depending on the current inputs configuration.

Process scripts can contain conditional statements by simply prefixing the script block with the keyword script:. By doing that the interpreter will evaluate all the following statements as a code block that must return the script string to be executed. It's much easier to use than to explain, for example:

```
seq_to_align = ...
mode = 'tcoffee'

process align {
    input:
    file seq_to_aln from sequences

    script:
    if( mode == 'tcoffee' )
        """
        t_coffee -in $seq_to_aln > out_file
        """

    else if( mode == 'mafft' )
        """
        mafft --anysymbol --parttree --quiet $seq_to_aln > out_file
        """

    else if( mode == 'clustalo' )
        """
        clustalo -i $seq_to_aln -o out_file
        """

    else
        error "Invalid alignment mode: ${mode}"
}
```
In the above example the process will execute the script fragment depending on the value of the mode parameter. By default it will execute the `tcoffee` command, changing the mode variable to `mafft` or `clustalo value, the other branches will be executed.

### Template

The Process script can be externalised by using template files which can be reused across different processes and tested independently from the overall pipeline execution.

A template is simply a shell script file that Nextflow is able to execute by using the template function as shown below:

```
process template_example {
    input:
    val STR from 'this', 'that'

    script:
    template 'my_script.sh'
}
```
Nextflow looks for the `my_script.sh` template file in the directory templates that must exist in the same folder where the Nextflow script file is located (any other location can be provided by using an absolute template path).

**NOTE:** When using [DSL2](https://www.nextflow.io/docs/latest/dsl2.html#dsl2-page), Nextflow also looks in the `templates` directory located in the same folder as module. See [module templates]().

### Shell

The `shell` block is a string statement that defines the *shell* command executed by the process to carry out its task. It is an alternative to the [Script](https://www.nextflow.io/docs/latest/process.html#process-script) definition with an important difference, it uses the exclamation mark `!` character as the variable placeholder for Nextflow variables in place of the usual dollar character.

In this way it is possible to use both Nextflow and Bash variables in the same piece of code without having to escape the latter and making process scripts more readable and easy to maintain. For example:

```
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process myTask {
  input:
  val str

  output:
  stdout

  shell:
  '''
  echo User $USER says !{str}
  '''
}

workflow {
  greetings_ch = Channel.of('Hello', 'Hola', 'Bonjour')
  myTask(greetings_ch).view()
}
```

In the above example the `$USER` variable is managed by the Bash interpreter, while `!{str}` is handled as a process input variable managed by Nextflow.

**NOTE**

* Shell script definitions require the use of single-quote `'` delimited strings. When using double-quote `"` delimited strings, dollar variables are interpreted as Nextflow variables as usual. See [String interpolation](https://www.nextflow.io/docs/latest/script.html#string-interpolation).
* Variables prefixed with `!` must always be enclosed in curly brackets, i.e. `!{str}` is a valid variable whereas `!str` is ignored.
* Shell scripts support the use of the file [Template](https://www.nextflow.io/docs/latest/process.html#process-template) mechanism. The same rules are applied to the variables defined in the script template.

### Native Execution

Nextflow processes can execute native code other than system scripts as shown in the previous paragraphs.

This means that instead of specifying the process command to be executed as a string script, you can define it by providing one or more language statements, as you would do in the rest of the pipeline script. Simply starting the script definition block with the `exec:` keyword, for example:

```
x = Channel.from( 'a', 'b', 'c')

process simpleSum {
    input:
    val x

    exec:
    println "Hello Mr. $x"
}
```

will display:

```
Hello Mr. b
Hello Mr. a
Hello Mr. c
```

## Stub

**WARNING:** This feature is experimental. It may change in future versions.

As of version 20.11.0-edge it's possible to define a command stub that replaces the actual process command, when the `-stub-run` or `-stub` command line option.

```
process INDEX {
  input:
    path transcriptome

  output:
    path 'index'

  script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """

  stub:
    """
    mkdir index
    touch index/seq.bin
    touch index/info.json
    touch index/refseq.bin
    """
}
```
This feature is meant to allow the fast prototyping and test of the workflow logic without using the real commands. The developer can use it to provide a dummy command which is expected to mimic the execution of the real one in a quicker manner. This can also be used as an alternative for the `dry-run` feature.

**TIP:** The `stub` block can be defined before or after the `script` block. When the pipeline is executed with the `-stub-run` option and a process's `stub` is not defined, the script block is executed.

## Inputs

Nextflow processes are isolated from each other but can communicate between themselves sending values through channels.

The `input` block defines from which channels the process expects to receive data. You can only define one input block at a time and it must contain one or more input declarations.

The input block follows the syntax shown below:

```
input:
  <input qualifier> <input name> [from <source channel>] [attributes]
```

An input definition starts with an input *qualifier* and the input name, followed by the keyword `from` and the actual channel over which inputs are received (this last part applies to DSL1). Finally some input optional attributes can be specified.

**TIP:** When the input name is the same as the channel name, the `from` part of the declaration can be omitted.

The input qualifier declares the type of data to be received. This information is used by Nextflow to apply the semantic rules associated to each qualifier and handle it properly depending on the target execution platform (grid, cloud, etc).

The qualifiers available are the ones listed in the following table:

| Qualifier   |  Semantic |
|-------------|-----------|
| val  | Lets you access the received input value by its name in the process script. |
| env  | Lets you use the received value to set an environment variable named as the specified input name. |
| file | Lets you handle the received value as a file, staging it properly in the execution context. |
| path | Lets you handle the received value as a path, staging the file properly in the execution context. |
| stdin| Lets you forward the received value to the process `stdin` special file. |
| tuple| Lets you handle a group of input values having one of the above qualifiers. |
| each | Lets you execute the process for each entry in the input collection. |

### Input of generic values

The `val` qualifier allows you to receive data of any type as input. It can be accessed in the process script by using the specified input name, as shown in the following example:

```
nextflow.enable.dsl=2

process basicExample {
  input:
  val x

  output:
  stdout

  shell:
  """
  echo process job !{x}
  """
}

workflow {
  num = Channel.of(1, 2, 3)
  basicExample(num).view()
}
```
In the above example the process is executed three times, each time a value is received from the channel `num` and used to process the script. Thus, it results in an output similar to the one shown below:

```
process job 3
process job 1
process job 2
```
**NOTE:** The *channel* guarantees that items are delivered in the same order as they were received - but - since the process is executed in a parallel manner, there is no guarantee that they are processed in the same order as they are received. In fact, in the above example, the value `3` is processed before the others.

### Input of files

The `file` qualifier allows the handling of file values in the process execution context. This means that Nextflow will stage it in the process execution directory, and it can be access in the script by using the name specified in the input declaration. For example:

```
nextflow.enable.dsl=2

process blastThemAll {
  input:
  path proteins

  shell:
  """
  blastp -query !{proteins} -db nr
  """
}

workflow {
  proteins = Channel.fromPath('/some/path/*.fa')
  blastThemAll(proteins)
}
```

In the above example all the files ending with the suffix `.fa` are sent over the channel proteins. Then, these files are received by the process which will execute a *BLAST* query on each of them.

There may be cases where your task needs to use a file whose name is fixed, it does not have to change along with the actual provided file. In this case you can specify its name by specifying the `name` attribute in the input file parameter declaration, as shown in the following example:

```
input:
    path query_file name 'query.fa' from proteins
```

Or alternatively using a shorter syntax:

```
input:
    path 'query.fa' from proteins
```

Using this, the previous example can be re-written as shown below:

```
proteins = Channel.fromPath( '/some/path/*.fa' )

process blastThemAll {
  input:
  path 'query.fa' from proteins

  "blastp -query query.fa -db nr"
}
```

What happens in this example is that each file, that the process receives, is staged with the name `query.fa` in a different execution context (i.e. the folder where the job is executed) and an independent process execution is launched.

**TIP:** This feature allows you to execute the process command multiple times without worrying about the file names changing. In other words, *Nextflow* helps you write pipeline tasks that are self-contained and decoupled from the execution environment. This is also the reason why you should avoid whenever possible using absolute or relative paths when referencing files in your pipeline processes.

### Multiple input files

A process can declare as input file a channel that emits a collection of values, instead of a simple value.

In this case, the script variable defined by the input file parameter will hold a list of files. You can use it as shown before, referring to all the files in the list, or by accessing a specific entry using the usual square brackets notation.

When a target file name is defined in the input parameter and a collection of files is received by the process, the file name will be appended by a numerical suffix representing its ordinal position in the list. For example:

```
#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2

process enumerateFasta {
  input:
  path "fasta_file"
  
  output:
  stdout
  
  shell:
  """
  echo fasta_file*
  """
}

workflow {
    fasta_ch = Channel.fromPath("fasta/*.fa").buffer(size:3)
    enumerateFasta(fasta_ch).view()
}
```
will output:

```
fasta_file1 fasta_file2 fasta_file3

fasta_file1 fasta_file2 fasta_file3

fasta_file1 fasta_file2 fasta_file3
...
```

The target input file name can contain the `*` and `?` wildcards, that can be used to control the name of staged files. The following table shows how the wildcards are replaced depending on the cardinality of the received input collection.

| Cardinality  | Name pattern  | Staged file names  |
|--------------|---------------|--------------------|
| any   | `*`  | named as the source file  |
| 1     | `file*.ext`  | `file.ext`  |
| 1     | `file?.ext`  | `file1.ext`  |
| 1     | `file??.ext`  | `file01.ext`  |
| many  | `file*.ext`  | `file1.ext`, `file2.ext`, `file3.ext`  |
| many  | `file?.ext`  | `file1.ext`, `file2.ext`, `file3.ext`  |
| many  | `file??.ext`  | `file01.ext`, `file02.ext`, `file03.ext`  |
| many  | `dir/*`  | named as the source file, created in `dir` subdirectory  |
| many  | `dir??/*`  | named as the source file, created in a progressively indexed subdirectory e.g. `dir01/`, `dir02/`, etc.  |
| many  | `dir*/*`  | as above  |

The following fragment shows how a wildcard can be used in the input file declaration:

```
fasta = Channel.fromPath( "/some/path/*.fa" ).buffer(size:3)

process blastThemAll {
    input:
    path 'seq?.fa' from fasta

    "cat seq1.fa seq2.fa seq3.fa"
}
```

**NOTE:** Rewriting input file names according to a named pattern is an extra feature and not at all obligatory. The normal file input constructs introduced in the [Input of files](https://www.nextflow.io/docs/latest/process.html#input-of-files) section are valid for collections of multiple files as well. To handle multiple input files preserving the original file names, use the `*` wildcard as name pattern or a variable identifier.

### Dynamic input file names

When the input file `name` is specified by using the name file clause or the short string notation, you are allowed to use other input values as variables in the file name string. For example:

```
process simpleCount {
  input:
  val x from species
  path "${x}.fa" from genomes

  """
  cat ${x}.fa | grep '>'
  """
}
```

In the above example, the input file name is set by using the current value of the `x` input value.

This allows the input files to be staged in the script working directory with a name that is coherent with the current execution context.

**TIP:** In most cases, you won’t need to use dynamic file names, because each process is executed in its own temporary directory, and input files are automatically staged into this directory by Nextflow. This guarantees that input files with the same name won’t overwrite each other.

### Input of type 'path'

The `path` input qualifier was introduced by Nextflow version 19.10.0 and it’s a drop-in replacement for the file qualifier, therefore it’s backward compatible with the syntax and the semantic for the input file described above.

The important difference between file and path qualifier is that the first expects the values received as input to be file objects. When inputs is a different type, it automatically coverts to a string and saves it to a temporary file. This can be useful in some uses cases, but it turned out to be tricky in most common cases.

The `path` qualifier instead interprets string values as the path location of the input file and automatically converts to a file object.

```
process foo {
  input:
    path x from '/some/data/file.txt'

  """
  your_command --in $x
  """
}
```

**NOTE:** The input value should represent an absolute path location, i.e. the string value must be prefixed with a `/` character or with a supported URI protocol (`file://`, `http://`, `s3://`, etc) and it cannot contain special characters (`\n`, etc).

The option `stageAs` allows you to control how the file should be named in the task work directory, providing a specific name or a name pattern as described in the [Multiple input files](https://www.nextflow.io/docs/latest/process.html#multiple-input-files) section:

```
process foo {
  input:
    path x, stageAs: 'file.txt' from '/some/data/file.txt'

  """
  your_command --in $x
  """
}
```

**TIP:** The `path` qualifier should be preferred over `file` to handle process input files when using Nextflow 19.10.0 or later.

### Input of type 'stdin'

The `stdin` input qualifier allows you the forwarding of the value received from a channel to the standard input of the command executed by the process. For example:

```
nextflow.enable.dsl=2

process printAll {
  input:
  stdin str

  output:
  stdout

  shell:
  """
  cat -
  """
}

workflow {
  greetings_ch = Channel.of('hello','hola', 'bonjour', 'ciao').map{ it + '\n'}
  printAll(greetings_ch).view()
}
```
it will output:

```
hola
bonjour
ciao
hello
```

### Input of type 'env'

The `env` qualifier allows you to define an environment variable in the process execution context based on the value received from the channel. For example:

```
nextflow.enable.dsl=2

process printEnv {
  input:
  env HELLO

  output:
  stdout

  '''
  echo $HELLO world!
  '''
}

workflow {
  str = Channel.of('hello', 'hola', 'bonjour', 'ciao')
  printEnv(str).view()
}
```

will output:

```
hello world!
ciao world!
bonjour world!
hola world!
```
### Input of type 'tuple'

The `tuple` qualifier allows you to group multiple parameters in a single parameter definition. It can be useful when a process receives, in input, tuples of values that need to be handled separately. Each element in the `tuple` is associated to a corresponding element with the tuple definition. For example:

```
values = Channel.of( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

process tupleExample {
    input:
    tuple val(x), path('latin.txt') from values

    """
    echo Processing $x
    cat - latin.txt > copy
    """
}
```

In the above example the `tuple` parameter is used to define the value `x` and the file `latin.txt`, which will receive a value from the same channel.

In the `tuple` declaration items can be defined by using the following qualifiers: `val`, `env`, `path` and `stdin`.

File names can be defined in *dynamic* manner as explained in the [Dynamic input file names](https://www.nextflow.io/docs/latest/process.html#dynamic-input-file-names) section.

### Input repeaters

The `each` qualifier allows you to repeat the execution of a process for each item in a collection, every time a new data is received. For example:

```
sequences = Channel.fromPath('*.fa')
methods = ['regular', 'expresso', 'psicoffee']

process alignSequences {
  input:
  path seq from sequences
  each mode from methods

  """
  t_coffee -in $seq -mode $mode > result
  """
}
```

In the above example every time a file of sequences is received as input by the process, it executes three tasks running a T-coffee alignment with a different value for the `mode` parameter. This is useful when you need to *repeat* the same task for a given set of parameters.

Input repeaters can be applied to files as well. For example:

```
sequences = Channel.fromPath('*.fa')
methods = ['regular', 'expresso']
libraries = [ file('PQ001.lib'), file('PQ002.lib'), file('PQ003.lib') ]

process alignSequences {
  input:
  path seq from sequences
  each mode from methods
  each path(lib) from libraries

  """
  t_coffee -in $seq -mode $mode -lib $lib > result
  """
}
```
**NOTE:** When multiple repeaters are declared, the process is executed for each combination of them.

In the latter example for any sequence input file emitted by the `sequences` channel are executed 6 alignments, 3 using the `regular` method against each library files, and other 3 by using the `expresso` method always against the same library files.

**TIP:** If you need to repeat the execution of a process over an n-tuple of elements instead of simple values or files, create a channel combining the input values as needed to trigger the process execution multiple times. Refer to the [combine](https://www.nextflow.io/docs/latest/operator.html#operator-combine), [cross](https://www.nextflow.io/docs/latest/operator.html#operator-cross) and [phase](https://www.nextflow.io/docs/latest/operator.html#operator-phase) operators for more details.

### Understand how multiple inputs work

A key feature of processes is the ability to handle inputs from multiple channels.

When two or more channels are declared as process inputs, the process stops until there’s a complete input configuration ie. it receives an input value from all the channels declared as input.

When this condition is verified, it consumes the input values coming from the respective channels, and spawns a task execution, then repeat the same logic until one or more channels have no more content.

This means channel values are consumed serially one after another and the first empty channel cause the process execution to stop even if there are other values in other channels.

For example:

```
nextflow.enable.dsl=2

process foo {
  input:
  val x
  val y

  output:
  stdout

  script:
  """
  echo $x and $y
  """
}

workflow {
  ch1 = Channel.of(1, 2)
  ch2 = Channel.of('a', 'b', 'c')
  foo(ch1, ch2).view()
}
```

The process `foo` is executed two times because the first input channel only provides two values and therefore the `c` element is discarded. It prints:

```
1 and a
2 and b
```

A different semantic is applied when using a *value channel* (a.k.a. *singleton channel*). This kind of channel is created by the [Channel.value](https://nextflow.io/docs/latest/channel.html#channel-value) factory method or implicitly when a process input specifies a simple value in the `from` clause. By definition, a value channel is bound to a single value and it can be read an unlimited number of times without consuming its content.

These properties make that when mixing a value channel with one or more (queue) channels, it does not affect the process termination because its content is applied repeatedly.

To better understand this behavior, compare the previous example with the following one:

```
nextflow.enable.dsl=2

process bar {
  input:
  val x
  val y

  output:
  stdout

  script:
  """
  echo $x and $y
  """
}

workflow {
  ch1 = Channel.value(1)
  ch2 = Channel.of('a', 'b', 'c')
  bar(ch1, ch2).view()
}
```

The above snippet executes the `bar` process three times because the first input is a value channel, therefore its content can be read as many times as needed. The process termination is determined by the content of the second channel. It prints:

```
1 and a
1 and b
1 and c
```

See also: [Channel types](https://nextflow.io/docs/latest/channel.html#channel-types).
## Outputs

The `output` declaration block allows you to define the channels used by the process to send out the results produced. You can only define one output block at a time and it must contain one or more output declarations.

The output block follows the syntax shown below:

```
output:
  <output qualifier> <output name> [into <target channel>[,channel,..]] [attribute [,..]]
```

Output definitions start by an output qualifier and the output name, followed by the keyword `into` and one or more channels over which outputs are sent (the latter part is not necessary in DSL2). Finally some optional attributes can be specified.

**TIP:** When the output name is the same as the channel name, the into part of the declaration can be omitted.

**NOTE:** If an output channel has not been previously declared in the pipeline script, it will be implicitly created by the output declaration itself.

The qualifiers that can be used in the output declaration block are the ones listed in the following table:

| Qualifier  | Semantic |
|------------|----------|
| val  | Sends variables with the name specified over the output channel. |
| file | Sends a file produced by the process with the name specified over the output channel. |
| path | Sends a file produced by the process with the name specified over the output channel (replaces file). |
| env  | Sends the variable defined in the process environment with the name specified over the output channel. |
| stdout | Sends the executed process `stdout` over the output channel. |
| tuple | Sends multiple values over the same output channel. |

### Output values

The `val` qualifier allows you to output a value defined in the script context. In a common usage scenario, this is a value which has been defined in the input declaration block, as shown in the following example:

```
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process foo {
    input:
    val x

    script:
    """
    echo $x > file."$x".txt
    """
}

workflow {
    methods = ['prot', 'dna', 'rna']
    method_ch = Channel.from(methods)
    method_ch.view{"Received: $it"}
    foo(method_ch)
}
```

Valid output values are value literals, input value identifiers, variables accessible in the process scope and value expressions. For example:

```
process foo {
  input:
  path fasta from 'dummy'

  output:
  val x into var_channel
  val 'BB11' into str_channel
  val "${fasta.baseName}.out" into exp_channel

  script:
  x = fasta.name
  """
  cat $x > file
  """
}
```

### Output files

The `file` qualifier allows you to output one or more files, produced by the process, over the specified channel. Now, the `path` qualifier should be used instead of `file`. For example:

```
// Enable DSL2
nextflow.enable.dsl=2

process randomNum {
  output:
  path 'result.txt'

  '''
  echo $RANDOM > result.txt
  '''
}

workflow {
    randomNum()
}
```

In the above example the process, when executed, creates a file named `result.txt` containing a random number.

### Dynamic output file names

When an output file name contains a `*` or `?` wildcard character it is interpreted as a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) path matcher. This allows you to *capture* multiple files into a list object and output them as a sole emission. For example:

```
// Enable DSL2
nextflow.enable.dsl=2

process splitLetters {
    output:
    path 'chunk_*'

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}

workflow {
  splitLetters()
    .flatMap()
    .subscribe{println "File: ${it.name} => ${it.text}"}
}
```

which prints:

```
File: chunk_aa => H
File: chunk_ab => o
File: chunk_ac => l
File: chunk_ad => a
```

**NOTE:** In the above example, the operator [flatMap](https://nextflow.io/docs/latest/operator.html#operator-flatmap) is used to transform the list of files emitted by the letters channel into a channel that emits each file object separately.

Some caveats on glob pattern behavior:

* Input files are not included (unless `includeInputs` is `true`)
* Directories are included, unless the `**` pattern is used to recurse through directories
  
**WARNING:** Although the input files matching a glob output declaration are not included in the resulting output channel, these files may still be transferred from the task scratch directory to the original task work directory. Therefore, to avoid unnecessary file copies, avoid using loose wildcards when defining output files, e.g. `file '*'`. Instead, use a prefix or a suffix to restrict the set of matching files to only the expected ones, e.g. file `'prefix_*.sorted.bam'`.

By default all the files matching the specified glob pattern are emitted by the channel as a sole (list) item. It is also possible to emit each file as a sole item by adding the `mode flatten` attribute in the output file declaration.

By using the `mode` attribute the previous example can be re-written as shown below:

```
// Enable DSL2
nextflow.enable.dsl=2

process splitLetters {
    output:
    path 'chunk_*'

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}

workflow {
  splitLetters()
    .flatten()
    .subscribe{println "File: ${it.name} => ${it.text}"}
}
```
Read more about glob syntax at the following link: [What is a glob?](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob)

### Output path

The `path` output qualifier was introduced by Nextflow version 19.10.0 and it’s a drop-in replacement for the `file` output qualifier, therefore it’s backward compatible with the syntax and the semantic for the input `file` described above.

The main advantage of `path` over the `file` qualifier is that it allows the specification of a number of options to fine-control the output files.

| Name  | Description |
|---|---|
| glob  | When `true` the specified name is interpreted as a glob pattern (default: `true`) |
| hidden  | When `true` hidden files are included in the matching output files (default: `false`)  |
| followLinks  | When `true` target files are returned in place of any matching symlink (default: `true`) |
| type  | Type of paths returned, either `file`, `dir` or `any` (default: `any`, or `file` if the specified file name pattern contains a `**` - double star - symbol)  |
| maxDepth  | Maximum number of directory levels to visit (default: *no limit*)  |
| includeInputs  | When `true` any input files matching an output file glob pattern are included. |

**WARNING:** The `file` qualifier interprets `:` as a path separator, therefore `file 'foo:bar'` captures two files named `foo` and `bar`. The `path` qualifier, on the other hand, does not, so the output definition `path 'foo:bar'` captures a single file named `foo:bar`.

**TIP:** The `path` qualifier should be preferred over `file` to handle process output files when using Nextflow 19.10.0 or later.

### Output 'stdout' special file

The `stdout` qualifier allows you to capture the stdout output of the executed process. For example:

```
// Enable DSL2
nextflow.enable.dsl=2

process sayHello {

  output:
  stdout

  script:
  """
  echo Hello world!
  """
}

workflow {
  sayHello().view()
}
```


### Output 'env'

The `env` qualifier allows you to capture a variable defined in the process execution environment and send it over the channel specified in the output parameter declaration:

```
// Enable DSL2
nextflow.enable.dsl=2

process myTask {
    output:
    env FOO

    script:
    '''
    FOO=$(ls -la)
    '''
}

workflow {
  myTask().view{"directory content: $it"}
}
```

### Output 'tuple' of values

The `tuple` qualifier allows you to send multiple values into a single channel. This feature is useful when you need to group together the results of multiple executions of the same process, as shown in the following example:

```
process blast {
  input:
  val species
  path query

  output:
  tuple val(species), path('result')

  script:
  """
  blast -db nr -query $query > result
  """
}

workflow {
  query_ch = Channel.fromPath('*.fa')
  species_ch = Channel.from('human', 'cow', 'horse')
  blast(species, query)
}
```

In the above example a `blast` task is executed for each pair of species and query that are received. When the task completes a new tuple containing the value for `species` and the file `result` is generated.

A `tuple` declaration can contain any combination of the following qualifiers, previously described: `val`, `path`, `env` and `stdout`.

File names can be defined in a dynamic manner as explained in the [Dynamic output file names](https://nextflow.io/docs/latest/process.html#process-dynoutname) section.

### Optional output

In most cases a process is expected to generate output that is added to the output channel. However, there are situations where it is valid for a process to not generate output. In these cases `optional true` may be added to the output declaration, which tells Nextflow not to fail the process if the declared output is not created.

```
output:
    path("output.txt") optional true
```

In this example, the process is normally expected to generate an `output.txt` file, but in the cases where the file is legitimately missing, the process does not fail. The output channel is only populated by those processes that do generate `output.txt`.

## When

The`when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value.

It is useful to enable/disable the process execution depending on the state of various inputs and parameters. For example:

```
// Enable DSL2
nextflow.enable.dsl=2

process find {
  input:
  path proteins
  val type

  when:
  proteins.name =~ /^BB11.*/ && type == 'nr'

  script:
  """
  blastp -query $proteins -db nr
  """
}

workflow {
  find()
}
```

## Directives
