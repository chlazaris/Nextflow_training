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

### groupBy (legacy)

**WARNING:** This operator is deprecated. Use the [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple) operator instead.

The `groupBy` operator collects the values emitted by the source channel and groups them together using a *mapping* function that associates each item with a key. When finished, it emits an associative array that maps each key to the set of items identified by that key.

For example:

```
Channel
    .from('hello', 'ciao', 'hola', 'hi', 'bonjour')
    .groupBy { String str -> str[0] }
    .view()
```

emits:

```
[ b:['bonjour'], c:['ciao'], h:['hello','hola','hi'] ]
```

The *mapping* function is an optional parameter. When omitted, the values are grouped according to these rules:

* Any value of type `Map` is associated with the value of its first entry, or null when the map itself is empty.
* Any value of type `Map.Entry` is associated with the value of its key attribute.
* Any value of type `Collection` or `Array` is associated with its first entry.
* For any other value, the value itself is used at the key.

### groupTuple

The `groupTuple` operator collects tuples (or lists) of values emitted by the source channel grouping together the elements that share the same key. Finally it emits a new tuple object for each distinct key collected.

In other words, the operator transforms a sequence of tuple like (*K, V, W, ..*) into a new channel emitting a sequence of (*K, list(V), list(W), ..*)

For example:

```
Channel
     .from( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
     .groupTuple()
     .view()
```

emits:

```
[1, [A, B, C]]
[2, [C, A]]
[3, [B, D]]
```

By default the first entry in the tuple is used as grouping key. A different key can be chosen by using the `by` parameter and specifying the index of the entry to be used as key (the index is zero-based). For example, grouping by the second value in each tuple:

```
Channel
     .from( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
     .groupTuple(by: 1)
     .view()
```

emits:

```
[[1, 2], A]
[[1, 3], B]
[[2, 1], C]
[[3], D]
```

Available parameters:

|Field|Description|
|-----|-----------|
|by|The index (zero based) of the element to be used as grouping key. A key composed by multiple elements can be defined specifying a list of indices e.g. by: `[0,2]`|
|sort|Defines the sorting criteria for the grouped items. See below for available sorting options.|
|size|The number of items the grouped list(s) has to contain. When the specified size is reached, the tuple is emitted.|
|remainder|When `false` incomplete tuples (i.e. with less than size grouped items) are discarded (default). When `true` incomplete tuples are emitted as the ending emission. Only valid when a `size` parameter is specified.|

Sorting options:

|Sort|Description|
|----|-----------|
|false|No sorting is applied (default).|
|true|Order the grouped items by the item natural ordering i.e. numerical for number, lexicographic for string, etc. See http://docs.oracle.com/javase/tutorial/collections/interfaces/order.html|
|hash|Order the grouped items by the hash number associated to each entry.|
|deep|Similar to the previous, but the hash number is created on actual entries content e.g. when the item is a file, the hash is created on the actual file content.|
|*custom*|A custom sorting criterion used to order the tuples element holding list of values. It can be specified by using either a [Closure](https://www.nextflow.io/docs/latest/script.html#script-closure) or a [Comparator](http://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html) object.
|

**TIP:** You should always specify the number of expected elements in each tuple using the size attribute to allow the `groupTuple` operator to stream the collected values as soon as possible. However, there are use cases in which each tuple has a different size depending on the grouping key. In this case use the built-in function `groupKey` that allows you to create a special grouping key object such that it’s possible to associate the group size for a given key.

### map

The `map` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. The function applied is called the *mapping* function and is expressed with a [closure](https://nextflow.io/docs/latest/script.html#script-closure) as shown in the example below:

```
Channel
    .from( 1, 2, 3, 4, 5 )
    .map { it * it }
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

emits:

```
1
4
9
16
25
Done
```

### reduce

The `reduce` operator applies a function of your choosing to every item emitted by a channel. Each time this function is invoked it takes two parameters: firstly the *i-th* emitted item and secondly the result of the previous invocation of the function itself. The result is passed on to the next function call, along with the *i+1* th item, until all the items are processed.

Finally, the `reduce` operator emits the result of the last invocation of your function as the sole output.

For example:

```
Channel
    .from( 1, 2, 3, 4, 5 )
    .reduce { a, b -> println "a: $a b: $b"; return a+b }
    .view { "result = $it" }
```

emits:

```
a: 1 b: 2
a: 3 b: 3
a: 6 b: 4
a: 10 b: 5
result = 15
```

**TIP:** A common use case for this operator is to use the first parameter as an *accumulator* the second parameter as the i-th item to be processed (THIS IS UNCLEAR).

Optionally you can specify a seed value in order to initialise the accumulator parameter as shown below:

```
myChannel.reduce( seedValue ) {  a, b -> ... }
```

### toList

The `toList` operator collects all the items emitted by a channel to a `List` object and emits the resulting collection as a single item. For example:

```
Channel
    .from( 1, 2, 3, 4 )
    .toList()
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

emits:

```
[1,2,3,4]
Done
```

See also: [collect](https://nextflow.io/docs/latest/operator.html#collect) operator.

### toSortedList

The `toSortedList` operator collects all the items emitted by a channel to a `List` object where they are sorted and emits the resulting collection as a single item. For example:

```
Channel
    .from( 3, 2, 1, 4 )
    .toSortedList()
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

emits:

```
[1,2,3,4]
Done
```

You may also pass a comparator closure as an argument to the `toSortedList` operator to customize the sorting criteria. For example, to sort by the second element of a tuple in descending order:

```
Channel
    .from( ["homer", 5], ["bart", 2], ["lisa", 10], ["marge", 3], ["maggie", 7])
    .toSortedList( { a, b -> b[1] <=> a[1] } )
    .view()
```

which emits:

```
[[lisa, 10], [maggie, 7], [homer, 5], [marge, 3], [bart, 2]]
```

See also: [collect](https://nextflow.io/docs/latest/operator.html#collect) operator.

### transpose

The `transpose` operator transforms a channel in such a way that the emitted items are the result of a transposition of all tuple elements in each item. For example:

```
Channel.from([
    ['a', ['p', 'q'], ['u','v'] ],
    ['b', ['s', 't'], ['x','y'] ]
    ])
    .transpose()
    .view()
```

emits:

```
[a, p, u]
[a, q, v]
[b, s, x]
[b, t, y]
```

Available parameters:

|Field|Description|
|----|----|
|by|The index (zero based) of the element to be transposed. Multiple elements can be defined specifying as list of indices e.g. `by: [0,2]|`
|remainder|When `false`, incomplete tuples are discarded (default). When `true`, incomplete tuples are emitted containing a `null` in place of a missing element.|

## Splitting operators

These operators are used to split items emitted by channels into chunks that can be processed by downstream operators or processes.

The available splitting operators are:

* [splitCsv](https://www.nextflow.io/docs/latest/operator.html#splitcsv)
* [splitFasta](https://www.nextflow.io/docs/latest/operator.html#splitfasta)
* [splitFastq](https://www.nextflow.io/docs/latest/operator.html#splitfastq)
* [splitText](https://www.nextflow.io/docs/latest/operator.html#splittext)

### splitCsv

The `splitCsv` operator allows you to parse text items emitted by a channel, that are formatted using the [CSV format](http://en.wikipedia.org/wiki/Comma-separated_values), and split them into records or group them into list of records with a specified length.

In the simplest case, just apply the `splitCsv` operator ro a channel emitting a CSV-formatted text file or text entry. For example:

```
Channel
    .from( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitCsv()
    .view(row -> "${row[0]} - ${row[1]} - ${row[2]}")
```

The above example shows how CSV text is parsed and is split into single rows. Values can be accessed by its column index in the row object.

When the CSV file begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its name, as shown in the following example:

```
Channel
    .from( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitCsv(header: true)
    .view { row -> "${row.alpha} - ${row.beta} - ${row.gamma}" }
```

which emits:

```
10 - 20 - 30
70 - 80 - 90
```

Alternatively, you can provide custom header names by specifying a the list of strings in the header parameter as shown below:

```
Channel
    .from( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitCsv(header: ['col1', 'col2', 'col3'], skip: 1 )
    .view { row -> "${row.col1} - ${row.col2} - ${row.col3}" }
```

Available parameters:

|Field|Description|
|-----|-----------|
|by|The number of rows in each *chunk*|
|sep|The character used to separate the values (default: `,`)|
|quote|Values may be quoted by single or double quote characters.|
|header|When `true` the first line is used as columns names. Alternatively it can be used to provide the list of columns names.|
|charset|Parse the content by using the specified charset e.g. `UTF-8`|
|strip|Removes leading and trailing blanks from values (default: `false`)|
|skip|Number of lines since the file beginning to ignore when parsing the CSV content.|
|limit|Limits the number of retrieved records for each file to the specified value.|
|decompress|When `true` decompress the content using the GZIP format before processing it (note: files whose name ends with .gz extension are decompressed automatically)|
|elem|The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element)|

### splitFasta

The `splitFasta` operator allows you to split the entries emitted by a channel, that are formatted using the [FASTA format](http://en.wikipedia.org/wiki/FASTA_format). It returns a channel which emits text item for each sequence in the received FASTA content.

The number of sequences in each text chunk produced by the `splitFasta` operator can be set by using the `by` parameter. The following example shows how to read a FASTA file and split it into chunks containing 10 sequences each:

```
Channel
     .fromPath('misc/sample.fa')
     .splitFasta( by: 10 )
     .view()
```

**WARNING:** Chunks are stored in memory by default. When splitting large files, specify the parameter `file: true` to save the chunks into files in order to avoid an `OutOfMemoryException`. See the parameter table below for details.

A second version of the `splitFasta` operator allows you to split a FASTA content into record objects, instead of text chunks. A record object contains a set of fields that let you access and manipulate the FASTA sequence information with ease.

In order to split a FASTA content into record objects, simply use the `record` parameter specifying the map of required the fields, as shown in the example below:

```
Channel
     .fromPath('misc/sample.fa')
     .splitFasta( record: [id: true, seqString: true ])
     .filter { record -> record.id =~ /^ENST0.*/ }
     .view { record -> record.seqString }
```

In this example, the file `misc/sample.fa` is split into records containing the `id` and the `seqString` fields (i.e. the sequence id and the sequence data). The following `filter` operator only keeps the sequences whose ID starts with the `ENST0` prefix, finally the sequence content is printed by using the `view` operator.

Available parameters:

|Field|Description|
|----|----|
|by|Defines the number of sequences in each *chunk* (default: `1`)|
|size|Defines the size in memory units of the expected chunks eg. 1.MB.|
|limit|Limits the number of retrieved sequences for each file to the specified value.|
|record|Parse each entry in the FASTA file as record objects (see following table for accepted values)|
|charset|Parse the content by using the specified charset e.g.` UTF-8`|
|compress|When `true` resulting file chunks are GZIP compressed. The `.gz` suffix is automatically added to chunk file names.|
|decompress|When `true`, decompress the content using the GZIP format before processing it (note: files whose name ends with .gz extension are decompressed automatically)|
|file|When true saves each split to a file. Use a string instead of true value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified folder.|
|elem|The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element)|

The following fields are available when using the `record` parameter:

|Field|Description|
|-----|-----------|
|id|The FASTA sequence identifier i.e. the word following the `>` symbol up to the first *blank* or *newline* character|
|header|The first line in a FASTA sequence without the `>` character|
|desc|The text in the FASTA header following the ID value|
|text|The complete FASTA sequence including the header|
|seqString|The sequence data as a single line string i.e. containing no *newline* characters|
|sequence|The sequence data as a multi-line string (always ending with a *newline* character)|
|width|Define the length of a single line when the `sequence` field is used, after that the sequence data continues on a new line.|

### splitFastq

### splitText

## Combining operators

## Forking operators

## Math operators

## Other operators