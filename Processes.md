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

### Conditional Scripts

### Template

### Shell

### Native Execution

## Stub

## Inputs

## Outputs

## When

## Directives
