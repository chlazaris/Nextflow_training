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

The `unique` operator allows you to remove duplicate items from a channel and only emit single items with no repetition.

For example:

```
Channel
    .from(1, 1, 1, 5, 7, 7, 7, 3, 3)
    .unique()
    .view()
```

```
1
5
7
3
```

You can also specify an optional [closure](https://www.nextflow.io/docs/latest/script.html#script-closure) that customizes the way it distinguishes between unique items. For example:

```
Channel
    .from(1, 2, 3, 5)
    .unique{ it % 2 }
    .view()
```

```
1
4
```

### distinct

The `distinct` operator allows you to remove *consecutive* duplicated items from a channel, so that each emitted item is different from the preceding one. For example:

```
Channel
    .from(1,1,2,2,2,3,1,1,2,2,3)
    .distinct()
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
2
3
1
2
3
Done
```

You can also specify an optional [closure](https://www.nextflow.io/docs/latest/script.html#script-closure) that customizes the way it distinguishes between distinct items. For example:

```
Channel
    .from(1,1,2,2,2,3,1,1,2,4,6)
    .distinct { it % 2 }
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
2
3
2
Done
```

### first

The `first` operator creates a channel that returns the first item emitted by the source channel, or eventually the first item that matches an optional condition. The condition can be specified by using a [regular expression](https://www.nextflow.io/docs/latest/script.html#script-regexp), a Java *class* type or any boolean *predicate*. For example:

```
// no condition is specified, emits the very first item: 1
Channel
    .from( 1, 2, 3 )
    .first()
    .view()

// emits the first String value: 'a'
Channel
    .from( 1, 2, 'a', 'b', 3 )
    .first( String )
    .view()

// emits the first item matching the regular expression: 'aa'
Channel
    .from( 'a', 'aa', 'aaa' )
    .first( ~/aa.*/ )
    .view()

// emits the first item for which the predicate evaluates to true: 4
Channel
    .from( 1,2,3,4,5 )
    .first { it > 3 }
    .view()
```

### randomSample

The `randomSample` operator allows you to create a channel emitting the specified number of items randomly taken from the channel to which is applied. For example:

```
Channel
      .from( 1..100 )
      .randomSample( 10 )
      .view()
```

The above snippet will print 10 numbers in the range from 1 to 100.

The operator supports a second parameter that allows you to set the initial seed for the random number generator. By setting it, the `randomSample` operator will always return the same pseudo-random sequence. For example:

```
Channel
      .from( 1..100 )
      .randomSample( 10, 234 )
      .view()
```

The above example will print 10 random numbers in the range between 1 and 100. At each run of the script, the same sequence will be returned.

### take

The `take` operator allows you to filter only the first *n* items emitted by a channel. For example:

```
Channel
    .from( 1,2,3,4,5,6 )
    .take( 3 )
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

```
1
2
3
Done
```

**TIP:** Specifying a size of `-1` causes the operator to take all values.

See also [until](https://www.nextflow.io/docs/latest/operator.html#until).

### last

The `last` operator creates a channel that only returns the last item emitted by the source channel. For example:

```
Channel
    .from(1, 2, 3, 4, 5, 6)
    .last()
    .view()
```

```
6
```

### until

The `until` operator creates a channel that returns the items emitted by the source channel and stop when the condition specified is verified (that last value is NOT included). For example:

```
Channel
    .from(3, 2, 1, 5, 1, 5)
    .until{ it == 5 }
    .view()
```

```
3
2
1
```

See also [take](https://www.nextflow.io/docs/latest/operator.html#take).

## Transforming operators

## Splitting operators

## Combining operators

## Forking operators

## Math operators

## Other operators