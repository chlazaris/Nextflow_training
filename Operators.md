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

Transforming operators are used to transform the items emitted by a channel to new values.

These operators are:

* [buffer](https://www.nextflow.io/docs/latest/operator.html#buffer)
* [collate](https://www.nextflow.io/docs/latest/operator.html#collate)
* [collect](https://www.nextflow.io/docs/latest/operator.html#collect)
* [flatten](https://www.nextflow.io/docs/latest/operator.html#flatten)
* [flatMap](https://www.nextflow.io/docs/latest/operator.html#flatmap)
* [groupBy](https://www.nextflow.io/docs/latest/operator.html#groupby)
* [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
* [map](https://www.nextflow.io/docs/latest/operator.html#map)
* [reduce](https://www.nextflow.io/docs/latest/operator.html#reduce)
* [toList](https://www.nextflow.io/docs/latest/operator.html#tolist)
* [toSortedList](https://www.nextflow.io/docs/latest/operator.html#tosortedlist)
* [transpose](https://www.nextflow.io/docs/latest/operator.html#transpose)

### buffer

The `buffer` operator gathers the items emitted by the source channel into subsets and emits these subsets separately.

There are a number of ways you can regulate how `buffer` gathers the items from the source channel into subsets:

* `buffer( closingCondition )`: starts to collect the items emitted by the channel into a subset until the closing condition is verified. After that the subset is emitted to the resulting channel and new items are gathered into a new subset. The process is repeated until the last value in the source channel is sent. The `closingCondition` can be specified either as a [regular expression](https://www.nextflow.io/docs/latest/script.html#script-regexp), a Java class, a literal value, or a *boolean predicate* that has to be satisfied. For example:

```
Channel
    .from( 1,2,3,1,2,3 )
    .buffer{ it==2 }
    .view()

// emitted values
[1,2]
[3,1,2]
```

* `buffer( openingCondition, closingCondition )`: starts to gather the items emitted by the channel as soon as one of the them verifies the *opening condition* and it continues until there is one item which verifies the *closing condition*. After that, the subset is emitted and it continues applying the described logic until the last channel item is emitted. Both conditions can be defined either as a [regular expression](https://www.nextflow.io/docs/latest/script.html#script-regexp), a literal value, a Java class, or a *boolean predicate* that need to be satisfied. For example:

```
Channel
    .from( 1,2,3,4,5,1,2,3,4,5,1,2 )
    .buffer( 2,4 )
    .view()

// emits bundles starting with '2' and ending with '4'
[2,3,4]
[2,3,4]
```

* `buffer( size: n )`: transforms the source channel in such a way that it emits tuples made up of `n` elements. An incomplete tuple is discarded. For example:

```
Channel
    .from( 1,2,3,1,2,3,1 )
    .buffer( size: 2 )
    .view()

// emitted values
[1, 2]
[3, 1]
[2, 3]
```

If you want to emit the last items in a tuple containing less than `n` elements, simply add the parameter `remainder` specifying `true`, for example:

```
Channel
    .from( 1,2,3,1,2,3,1 )
    .buffer( size: 2, remainder: true )
    .view()

// emitted values
[1, 2]
[3, 1]
[2, 3]
[1]
```

* `buffer( size: n, skip: m )`: as in the previous example, it emits tuples containing `n` elements, but skips `m` values before starting to collect the values for the next tuple (including the first emission). For example:

```
Channel
    .from( 1,2,3,4,5,1,2,3,4,5,1,2 )
    .buffer( size:3, skip:2 )
    .view()

// emitted values
[3, 4, 5]
[3, 4, 5]
```

If you want to emit the remaining items in a tuple containing less than `n` elements, simply add the parameter `remainder` specifying `true`, as shown in the previous example.

See also: [collate](https://www.nextflow.io/docs/latest/operator.html#collate) operator.

### collate

The `collate` operator transforms a channel in such a way that the emitted values are grouped in tuples containing n items. For example:

```
Channel
    .from(1,2,3,1,2,3,1)
    .collate( 3 )
    .view()
```

emits:

```
[1, 2, 3]
[1, 2, 3]
[1]
```

As shown in the above example the last tuple may be incomplete e.g. contain fewer elements than the specified size. If you want to avoid this, specify `false` as the second parameter. For example:

```
Channel
    .from(1,2,3,1,2,3,1)
    .collate( 3, false )
    .view()
```

emits:

```
[1, 2, 3]
[1, 2, 3]
```

A second version of the `collate` operator allows you to specify, after the *size*, the *step* by which elements are collected in tuples. For example:

```
Channel
    .from(1,2,3,4)
    .collate( 3, 1 )
    .view()
```

emits:

```
[1, 2, 3]
[2, 3, 4]
[3, 4]
[4]
```

As before, if you don’t want to emit the last items which do not complete a tuple, specify `false` as the third parameter.

See also: [buffer](https://www.nextflow.io/docs/latest/operator.html#buffer) operator.

### collect

The `collect` operator collects all the items emitted by a channel to a `List` and returns the resulting object as a sole emission. For example:

```
Channel
    .from( 1,2,3,4 )
    .collect()
    .view()

# emits
[1,2,3,4]
```

An optional [closure]() can be specified to transform each item before adding it to the resulting list. For example:

```
Channel
    .from( 'hello','ciao','bonjour' )
    .collect{ it.length() }
    .view()

# emits
[5,4,7]
```

See also: [toList](https://www.nextflow.io/docs/latest/operator.html#tolist) and [toSortedList](https://www.nextflow.io/docs/latest/operator.html#tosortedlist) operator.

### flatten

The `flatten` operator transforms a channel in such a way that every item of type `Collection` or `Array` is flattened so that each single entry is emitted separately by the resulting channel. For example:

```
Channel
    .from( [1,[2,3]], 4, [5,[6]] )
    .flatten()
    .view()
```

emits: 

```
1
2
3
4
5
6
```

See also: `flatMap` operator
### flatMap

The `flatMap` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. Whereas the mapping function returns a list of items, this list is flattened so that each single item is emitted on its own.

For example:

```
// create a channel of numbers
numbers = Channel.from( 1, 2, 3 )

// map each number to a tuple (array), which items are emitted separately
results = numbers.flatMap { n -> [ n*2, n*3 ] }

// print the final results
results.subscribe onNext: { println it }, onComplete: { println 'Done' }
```

emits:

```
2
3
4
6
6
9
Done
```

Associative arrays are handled in the same way, so that each array entry is emitted as a single key-value item. For example:

```
Channel
    .from ( 1, 2, 3 )
    .flatMap { it -> [ number: it, square: it*it ] }
    .view { it.key + ': ' + it.value }
```

emits:

```
number: 1
square: 1
number: 2
square: 4
number: 3
square: 9
```
## Splitting operators

## Combining operators

## Forking operators

## Math operators

## Other operators