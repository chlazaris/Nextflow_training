# Operators

Nextflow *operators* are methods that allow you to connect channels to each other or to transform values emitted by a channel applying some user provided rules.

Operators can be separated into seven groups:

* [Filtering operators](https://www.nextflow.io/docs/latest/operator.html#filtering-operators)
* [Transforming operators](https://www.nextflow.io/docs/latest/operator.html#transforming-operators)
* [Splitting operators](https://www.nextflow.io/docs/latest/operator.html#splitting-operators)
* [Combining operators](https://www.nextflow.io/docs/latest/operator.html#combining-operators)
* [Forking operators](https://www.nextflow.io/docs/latest/operator.html#forking-operators)
* [Math operators](https://www.nextflow.io/docs/latest/operator.html#maths-operators)
* [Other operators](https://www.nextflow.io/docs/latest/operator.html#other-operators)

**NOTE:** The operators [set](https://www.nextflow.io/docs/latest/operator.html#operator-set) and `subscribe` are *final* operators and therefore, if used, they must be the last operator in a chain of combined operators.

## Filtering operators

Given a channel, filtering operators allow you to select only the items that comply with a given rule.

The available filtering operators are:

* [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
* [unique](https://www.nextflow.io/docs/latest/operator.html#unique)
* [distinct](https://www.nextflow.io/docs/latest/operator.html#distinct)
* [first](https://www.nextflow.io/docs/latest/operator.html#first)
* [randomSample](https://www.nextflow.io/docs/latest/operator.html#randomsample)
* [take](https://www.nextflow.io/docs/latest/operator.html#take)
* [last](https://www.nextflow.io/docs/latest/operator.html#last)
* [until](https://www.nextflow.io/docs/latest/operator.html#until)

### filter

The `filter` operator allows you to get only the items emitted by a channel that satisfy a condition, discarding all the others. The filtering condition can be specified by using either a [regular expression](https://www.nextflow.io/docs/latest/script.html#script-regexp), a literal value, a type *qualifier* (i.e. a Java class) or any boolean *predicate*

The following example shows how to filter a channel by using a regular expression that returns only strings that begin with `a`:

```
Channel
    .from('a', 'b', 'aa', 'bc', 3, 4.5)
    .filter(~/^a.*/)
    .view()
```
```
a
aa
```

The following example shows how to filter a channel by specifying the type qualifier `Number` so that only numbers are returned:

```
Channel
    .from('a', 'b', 'aa', 'bc', 3, 4.5)
    .filter( Number )
    .view()
```

```
3
4.5
```

Finally, a filtering condition can be defined by using any a boolean *predicate*. A predicate is expressed by a [closure](https://www.nextflow.io/docs/latest/script.html#script-closure) returning a boolean value. For example the following fragment shows how filter a channel emitting numbers so that the *odd* values are returned:

```
Channel
    .from(1, 2, 3, 4, 5)
    .filter{it % 2 == 1}
    .view()
```

```
1
3
5
```

**TIP:** In the above example the filter condition is wrapped in curly brackets, instead of parentheses, because it specifies a [closure](https://www.nextflow.io/docs/latest/script.html#script-closure) as the operator’s argument. In reality it is just syntactic sugar for
`filter({it % 2 == 1})`

### unique

### distinct

### first

### randomSample

### take

### last

### until

## Transforming operators

## Splitting operators

## Combining operators

## Forking operators

## Math operators

## Other operators