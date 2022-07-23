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

## Outputs

## When

## Directives
