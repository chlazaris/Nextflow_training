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

The `splitFastq` operator allows you to split the entries emitted by a channel, that are formatted using the [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format). It returns a channel which emits a text chunk for each sequence in the received item.

The number of sequences in each text chunk produced by the `splitFastq` operator is defined by the parameter `by`. The following example shows you how to read a FASTQ file and split it into chunks containing 10 sequences each:

```
Channel
    .fromPath('misc/sample.fastq')
    .splitFastq( by: 10 )
    .view()
```

**WARNING:** Chunks are stored in memory by default. When splitting large files, specify the parameter `file: true` to save the chunks into files in order to avoid an `OutOfMemoryException`. See the parameter table below for details.

A second version of the `splitFastq` operator allows you to split a FASTQ formatted content into record objects, instead of text chunks. A record object contains a set of fields that let you access and manipulate the FASTQ sequence data with ease.

In order to split FASTQ sequences into record objects simply use the `record` parameter specifying the map of the required fields, or just specify `record: true` as in the example shown below:

```
Channel
    .fromPath('misc/sample.fastq')
    .splitFastq( record: true )
    .view { record -> record.readHeader }
```

Finally the `splitFastq` operator is able to split paired-end read pair FASTQ files. It must be applied to a channel which emits tuples containing at least two elements that are the files to be splitted. For example:

```
Channel
    .fromFilePairs('/my/data/SRR*_{1,2}.fastq', flat: true)
    .splitFastq(by: 100_000, pe: true, file: true)
    .view()
```

**NOTE:** The `fromFilePairs` requires the `flat: true` option in order to emit the file pairs as separate elements in the produced tuples.

**NOTE:** This operator assumes that the order of the paired-end reads correspond with each other and both files contain the same number of reads.

Available parameters:

|Field|Description|
|----|----|
|by|Defines the number of reads in each *chunk* (default: `1`)|
|pe|When `true` splits paired-end read files, therefore items emitted by the source channel must be tuples in which at least two elements are the read-pair files to be splitted.|
|limit|Limits the number of retrieved *reads* for each file to the specified value.|
|record|Parse each entry in the FASTQ file as record objects (see following table for accepted values)|
|charset|Parse the content by using the specified charset e.g.` UTF-8`|
|compress|When `true` resulting file chunks are GZIP compressed. The `.gz` suffix is automatically added to chunk file names.|
|decompress|When `true`, decompress the content using the GZIP format before processing it (note: files whose name ends with .gz extension are decompressed automatically)|
|file|When true saves each split to a file. Use a string instead of true value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in order to save the split files into the specified folder.|
|elem|The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element)|

The following fields are available when using the `record` parameter:

|Field|Description|
|-----|-----------|
|readHeader|Sequence header (without the `@` prefix)|
|readString|The raw sequence data|
|qualityHeader|Base quality header (it may be empty)|
|qualityString|Quality values for the sequence|

### splitText

The `splitText` operator allows you to split multi-line strings or text file items, emitted by a source channel into chunks containing n lines, which will be emitted by the resulting channel.

For example:

```
Channel
    .fromPath('/some/path/*.txt')
    .splitText()
    .view()
```

It splits the content of the files with suffix `.txt`, and prints it line by line.

By default, the `splitText` operator splits each item into chunks of one line. You can define the number of lines in each chunk by using the parameter `by`, as shown in the following example:

```
Channel
    .fromPath('/some/path/*.txt')
    .splitText( by: 10 )
    .subscribe {
        print it;
        print "--- end of the chunk ---\n"
    }
```

An optional [closure](https://www.nextflow.io/docs/latest/script.html#script-closure) can be specified in order to transform the text chunks produced by the operator. The following example shows how to split text files into chunks of 10 lines and transform them to capital letters:

```
Channel
    .fromPath('/some/path/*.txt')
    .splitText( by: 10 ) { it.toUpperCase() }
    .view()
```

**NOTE:** Text chunks returned by the operator `splitText` are always terminated by a `\n` newline character.

Available parameters:

|Field|Description|
|-----|-----|
|by|Defines the number of lines in each chunk (default: `1`).|
|limit|Limits the number of retrieved lines for each file to the specified value.|
|charset|Parse the content by using the specified charset e.g. `UTF-8.`|
|compress|When `true` resulting file chunks are GZIP compressed. The `.gz` suffix is automatically added to chunk file names.|
|decompress|When true, decompress the content using the GZIP format before processing it (note: files whose name ends with `.gz` extension are decompressed automatically).|
|file|When `true` saves each split to a file. Use a string instead of `true` value to create split files with a specific name (split index number is automatically added). Finally, set this attribute to an existing directory, in oder to save the split files into the specified folder.|
|elem|The index of the element to split when the operator is applied to a channel emitting list/tuple objects (default: first file object or first element).|
|keepHeader|Parses the first line as header and prepends it to each emitted chunk.|

## Combining operators

The combining operators are:

* [cross](https://www.nextflow.io/docs/latest/operator.html#cross)
* [collectFile](https://www.nextflow.io/docs/latest/operator.html#collectfile)
* [combine](https://www.nextflow.io/docs/latest/operator.html#combine)
* [concat](https://www.nextflow.io/docs/latest/operator.html#concat)
* [join](https://www.nextflow.io/docs/latest/operator.html#join)
* [merge](https://www.nextflow.io/docs/latest/operator.html#merge)
* [mix](https://www.nextflow.io/docs/latest/operator.html#mix)
* [phase](https://www.nextflow.io/docs/latest/operator.html#phase)
* [spread](https://www.nextflow.io/docs/latest/operator.html#spread)
* [tap](https://www.nextflow.io/docs/latest/operator.html#tap)

### cross

The `cross` operator allows you to combine the items of two channels in such a way that the items of the source channel are emitted along with the items emitted by the target channel for which they have a matching key.

The key is defined, by default, as the first entry in an array, a list or map object, or the value itself for any other data type. For example:

```
source = Channel.from( [1, 'alpha'], [2, 'beta'] )
target = Channel.from( [1, 'x'], [1, 'y'], [1, 'z'], [2,'p'], [2,'q'], [2,'t'] )

source.cross(target).view()
```

emits:

```
[ [1, alpha], [1, x] ]
[ [1, alpha], [1, y] ]
[ [1, alpha], [1, z] ]
[ [2, beta],  [2, p] ]
[ [2, beta],  [2, q] ]
[ [2, beta],  [2, t] ]
```

The above example shows how the items emitted by the source channels are associated to the ones emitted by the target channel (on the right) having the same key.

There are two important caveats when using the `cross` operator:

1. The operator is not *commutative*, i.e. the result of `a.cross(b)` is different from `b.cross(a)`
2. The source channel emits items for which there’s no key repetition i.e. the emitted items have an unique key identifier.

Optionally, a mapping function can be specified in order to provide a custom rule to associate an item to a key, in a similar manner as shown for the [phase](https://www.nextflow.io/docs/latest/operator.html#phase) operator.

### collectFile

The `collectFile` operator allows you to gather the items emitted by a channel and save them to one or more files. The operator returns a new channel that emits the collected file(s).

In the simplest case, just specify the name of a file where the entries have to be stored. For example:

```
Channel
    .from('alpha', 'beta', 'gamma')
    .collectFile(name: 'sample.txt', newLine: true)
    .subscribe {
        println "Entries are saved to file: $it"
        println "File content is: ${it.text}"
    }
```

A second version of the `collectFile` operator allows you to gather the items emitted by a channel and group them together into files whose name can be defined by a dynamic criteria. The grouping criteria is specified by a [closure](https://www.nextflow.io/docs/latest/script.html#script-closure) that must return a pair in which the first element defines the file name for the group and the second element the actual value to be appended to that file. For example:

```
Channel
    .from('Hola', 'Ciao', 'Hello', 'Bonjour', 'Halo')
    .collectFile() { item ->
        [ "${item[0]}.txt", item + '\n' ]
    }
    .subscribe {
        println "File ${it.name} contains:"
        println it.text
    }
```

emits:

```
File 'B.txt' contains:
Bonjour

File 'C.txt' contains:
Ciao

File 'H.txt' contains:
Halo
Hola
Hello
```

**TIP:** When the items emitted by the source channel are files, the grouping criteria can be omitted. In this case the items content will be grouped into file(s) having the same name as the source items.

The following parameters can be used with the `collectFile` operator:

|Name|Description|
|----|-----------|
|`cache`|Controls the caching ability of the `collectFile` operator when using the *resume* feature. It follows the same semantic of the [cache](https://www.nextflow.io/docs/latest/process.html#process-cache) directive (default: `true`).|
|`keepHeader`|Prepend the resulting file with the header fetched in the first collected file. The header size (ie. lines) can be specified by using the `skip` parameter (default: `false`), to determine how many lines to remove from all collected files except for the first (where no lines will be removed).|
|`name`|Name of the file where all received values are stored.|
|`newLine`|Appends a `newline` character automatically after each entry (default: `false`).|
|`seed`|A value or a map of values used to initialise the files content.|
|`skip`|Skip the first n lines eg. `skip: 1`.|
|`sort`|Defines sorting criteria of content in resulting file(s). See below for sorting options.|
|`storeDir`|Folder where the resulting file(s) are be stored.|
|`tempDir`|Folder where temporary files, used by the collecting process, are stored.|

**NOTE:** The file content is sorted in such a way that it does not depend on the order in which entries were added to it, which guarantees that it is consistent (i.e. does not change) across different executions with the same data.

The ordering of file’s content can be defined by using the `sort` parameter. The following criteria can be specified:

|Sort|Description|
|----|-----------|
|`false`|Disable content sorting. Entries are appended as they are produced.|
|`true`|Order the content by the entries natural ordering i.e. numerical for number, lexicographic for string, etc. See: http://docs.oracle.com/javase/tutorial/collections/interfaces/order.html|
|`index`|Order the content by the incremental index number assigned to each entry while they are collected.|
|`hash`|Order the content by the hash number associated to each entry (default)|
|`deep`|Similar to the previous, but the hash number is created on actual entries content e.g. when the entry is a file the hash is created on the actual file content.|
|*custom*|A custom sorting criteria can be specified by using either a [Closure](https://www.nextflow.io/docs/latest/script.html#script-closure) or a [Comparator](http://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html) object.|

For example the following snippet shows how sort the content of the result file alphabetically:

```
Channel
    .from('Z'..'A')
    .collectFile(name:'result', sort: true, newLine: true)
    .view { it.text }
```

emits:

```
A
B
C
:
Z
```

The following example shows how use a closure to collect and sort all sequences in a FASTA file from shortest to longest:

```
Channel
    .fromPath('/data/sequences.fa')
    .splitFasta( record: [id: true, sequence: true] )
    .collectFile( name:'result.fa', sort: { it.size() } )  {
        it.sequence
    }
    .view { it.text }
```

**WARNING:** The `collectFile` operator needs to store files in a temporary folder that is automatically deleted on workflow completion. For performance reasons this folder is located in the machine’s local storage, and it will require as much free space as the data that is being collected. Optionally, a different temporary data folder can be specified by using the `tempDir` parameter.

### combine

The combine `operator` combines (cartesian product) the items emitted by two channels or by a channel and a `Collection` object (as right operand). For example:

```
numbers = Channel.from(1,2,3)
words = Channel.from('hello', 'ciao')
numbers
    .combine(words)
    .view()
```

emits:

```
[1, hello]
[2, hello]
[3, hello]
[1, ciao]
[2, ciao]
[3, ciao]
```

A second version of the combine operator allows you to combine between them those items that share a common matching key. The index of the key element is specified by using the by parameter (the index is zero-based, multiple indexes can be specified with list a integers). For example:

```
left = Channel.from(['A',1], ['B',2], ['A',3])
right = Channel.from(['B','x'], ['B','y'], ['A','z'], ['A', 'w'])

left
    .combine(right, by: 0)
    .view()
```

emits:

```
[A, 1, z]
[A, 3, z]
[A, 1, w]
[A, 3, w]
[B, 2, x]
[B, 2, y]
```

See also [join](https://www.nextflow.io/docs/latest/operator.html#join), [cross](https://www.nextflow.io/docs/latest/operator.html#cross), [spead](https://www.nextflow.io/docs/latest/operator.html#spread), and [phase](https://www.nextflow.io/docs/latest/operator.html#phase)

### concat

The `concat` operator allows you to concatenate the items emitted by two or more channels to a new channel, in such a way that the items emitted by the resulting channel are in same order as they were when specified as operator arguments.

In other words it guarantees that given any *n* channels, the concatenation channel emits the items proceeding from the channel *i+1 th* only after all the items proceeding from the channel *i th* were emitted.

For example:

```
a = Channel.from('a','b','c')
b = Channel.from(1,2,3)
c = Channel.from('p','q')

c.concat( b, a ).view()
```

emits:

```
p
q
1
2
3
a
b
c
```

### join

The `join` operator creates a channel that joins together the items emitted by two channels for which exists a matching key. The key is defined, by default, as the first element in each item emitted.

For example:

```
left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right= Channel.from(['Z', 6], ['Y', 5], ['X', 4])
left.join(right).view()
```

emits:

```
[Z, 3, 6]
[Y, 2, 5]
[X, 1, 4]
```

The *index* of a different matching element can be specified by using the `by` parameter.

The `join` operator can emit all the pairs that are incomplete, i.e. the items for which a matching element is missing, by specifying the optional parameter `remainder` as shown below:

```
left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right= Channel.from(['Z', 6], ['Y', 5], ['X', 4])
left.join(right, remainder: true).view()
```

emits:

```
[Y, 2, 5]
[Z, 3, 6]
[X, 1, 4]
[P, 7, null]
```

The following parameters can be used with the join operator:

|Name|Description|
|----|-----------|
|by|The index (zero based) of the element to be used as grouping key. A key composed by multiple elements can be defined specifying a list of indices e.g. `by: [0,2]`|
|remainder|When `false` incomplete tuples (i.e. with less than size grouped items) are discarded (default). When `true` incomplete tuples are emitted as the ending emission.|
|failOnDuplicate|An error is reported when the same key is found more than once.|
|failOnMismatch|An error is reported when a channel emits a value for which there isn’t a corresponding element in the joining channel. This option cannot be used with `remainder`.|

### merge

The `merge` operator lets you join items emitted by two (or more) channels into a new channel.

For example the following code merges two channels together, one which emits a series of odd integers and the other which emits a series of even integers:

```
odds  = Channel.from([1, 3, 5, 7, 9]);
evens = Channel.from([2, 4, 6]);

odds
    .merge( evens )
    .view()
```

emits:

```
[1, 2]
[3, 4]
[5, 6]
```

An option closure can be provide to customize the items emitted by the resulting merged channel. For example:

```
odds  = Channel.from([1, 3, 5, 7, 9]);
evens = Channel.from([2, 4, 6]);

odds
    .merge( evens ) { a, b -> tuple(b*b, a) }
    .view()
```

emits:

```
[4, 1]
[16, 3]
[36, 5]
```

**NOTE** If using DSL2, you will receive a message that `merge` is deprecated and will be removed in a future release.

**DANGER:** When this operator is used to *merge* the outputs of two processes, keep in mind that the resulting merged channel will have non-deterministic behavior and may cause your pipeline execution to not resume properly. Because each process is executed in parallel and produces its outputs independently, there is no guarantee that they will be executed in the same order. Therefore the content of the resulting merged channel may have a different order on each run and may cause the resume to not work properly. For a better alternative use the [join](https://www.nextflow.io/docs/latest/operator.html#join) operator instead.

### mix

The `mix` operator combines the items emitted by two (or more) channels into a single channel.

For example:

```
c1 = Channel.from( 1,2,3 )
c2 = Channel.from( 'a','b' )
c3 = Channel.from( 'z' )

c1.mix(c2,c3)
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

emits:

```
1
2
3
'a'
'b'
'z'
```

**NOTE:** The items emitted by the resulting mixed channel may appear in any order, regardless of which source channel they came from. Thus, the following example could also be a possible result of the above example:

```
'z'
1
'a'
2
'b'
3
```

To make sure that the order is retained, use [concat](https://www.nextflow.io/docs/latest/operator.html#concat).

### phase

**WARNING:** This operator is deprecated. Use the [join][(https://www.nextflow.io/docs/latest/operator.html#join) operator instead.

The `phase` operator creates a channel that synchronizes the values emitted by two other channels, in such a way that it emits pairs of items that have a matching key.

The key is defined, by default, as the first entry in an array, a list or map object, or the value itself for any other data type.

For example:

```
ch1 = Channel.from( 1,2,3 )
ch2 = Channel.from( 1,0,0,2,7,8,9,3 )
ch1 .phase(ch2) .view()
```

emits:

```
[1,1]
[2,2]
[3,3]
```

Optionally, a mapping function can be specified in order to provide a custom rule to associate an item to a key, as shown in the following example:

```
ch1 = Channel.from( [sequence: 'aaaaaa', id: 1], [sequence: 'bbbbbb', id: 2] )
ch2 = Channel.from( [val: 'zzzz', id: 3], [val: 'xxxxx', id: 1], [val: 'yyyyy', id: 2])
ch1 .phase(ch2) { it -> it.id } .view()
```

emits:

```
[[sequence:aaaaaa, id:1], [val:xxxxx, id:1]]
[[sequence:bbbbbb, id:2], [val:yyyyy, id:2]]
```

Finally, the `phase` operator can emit all the pairs that are incomplete, i.e. the items for which a matching element is missing, by specifying the optional parameter remainder as shown below:

```
ch1 = Channel.from( 1,0,0,2,5,3 )
ch2 = Channel.from( 1,2,3,4 )
ch1 .phase(ch2, remainder: true) .view()
```

emits:

```
[1, 1]
[2, 2]
[3, 3]
[0, null]
[0, null]
[5, null]
[null, 4]
```

See also the [join](https://www.nextflow.io/docs/latest/operator.html#join) operator.

### spread

**WARNING:** This operator is deprecated. Use the [combine](https://www.nextflow.io/docs/latest/operator.html#combine) operator instead.

The `spread` operator combines the items emitted by the source channel with all the values in an array or a `Collection` object specified as the operator argument. For example:

```
Channel
    .from(1,2,3)
    .spread(['a','b'])
    .subscribe onNext: { println it }, onComplete: { println 'Done' }
```

emits:

```
[1, 'a']
[1, 'b']
[2, 'a']
[2, 'b']
[3, 'a']
[3, 'b']
Done
```

### tap

The `tap` operator combines the functions of [into](https://www.nextflow.io/docs/latest/operator.html#into) and [separate](https://www.nextflow.io/docs/latest/operator.html#separate) operators in such a way that it connects two channels, copying the values from the source into the *tapped* channel. At the same time it splits the source channel into a newly created channel that is returned by the operator itself.

The `tap` operator can be useful in certain scenarios where you may be required to concatenate multiple operations, as in the following example:

```
log1 = Channel.create()
log2 = Channel.create()

Channel
    .of ( 'a', 'b', 'c' )
    .tap ( log1 )
    .map { it * 2 }
    .tap ( log2 )
    .map { it.toUpperCase() }
    .view { "Result: $it" }

log1.view { "Log 1: $it" }
log2.view { "Log 2: $it" }
```

which emits:

```
Result: AA
Result: BB
Result: CC

Log 1: a
Log 1: b
Log 1: c

Log 2: aa
Log 2: bb
Log 2: cc
```

The `tap` operator also allows the target channel to be specified by using a closure. The advantage of this syntax is that you won’t need to previously create the target channel, because it is created implicitly by the operator itself.

Using the closure syntax the above example can be rewritten as shown below:

```
Channel
    .of ( 'a', 'b', 'c' )
    .tap { log1 }
    .map { it * 2 }
    .tap { log2 }
    .map { it.toUpperCase() }
    .view { "Result: $it" }

log1.view { "Log 1: $it" }
log2.view { "Log 2: $it" }
```

See also [into](https://www.nextflow.io/docs/latest/operator.html#into) and [separate](https://www.nextflow.io/docs/latest/operator.html#separate) operators.

## Forking operators

The forking operators are:

* [branch](https://www.nextflow.io/docs/latest/operator.html#branch)
* [choice](https://www.nextflow.io/docs/latest/operator.html#choice)
* [multiMap](https://www.nextflow.io/docs/latest/operator.html#multimap)
* [into](https://www.nextflow.io/docs/latest/operator.html#into)
* [separate](https://www.nextflow.io/docs/latest/operator.html#separate)
* [tap](https://www.nextflow.io/docs/latest/operator.html#tap)

### branch

**NOTE:** Requires Nextflow version `19.08.0-edge` or later.

The `branch` operator allows you to forward the items emitted by a source channel to one or more output channels, *choosing* one out of them at a time.

The selection criteria are defined by specifying a [closure](https://www.nextflow.io/docs/latest/script.html#script-closure) that provides one or more boolean expressions, each of which is identified by a unique label. On the first expression that evaluates to a true value, the current item is bound to a named channel as the label identifier. For example:

```
Channel
    .from(1,2,3,40,50)
    .branch {
        small: it < 10
        large: it > 10
    }
    .set { result }

 result.small.view { "$it is small" }
 result.large.view { "$it is large" }
```

emits:

```
1 is small
2 is small
3 is small
40 is large
50 is large
```

**NOTE:** The above *small* and *large* strings may be printed in any order due to the asynchronous execution of the `view` operator.

A default fallback condition can be specified using `true` as the last branch condition:

```
Channel
    .from(1,2,3,40,50)
    .branch {
        small: it < 10
        large: it < 50
        other: true
    }
```

The value returned by each branch condition can be customized by specifying an optional expression statement(s) just after the condition expression. For example:

```
Channel
    .from(1,2,3,40,50)
    .branch {
        foo: it < 10
            return it+2

        bar: it < 50
            return it-2

        other: true
            return 0
    }
```

**TIP:** When the `return` keyword is omitted, the value of the last statement is implicitly returned.

To create a branch criteria as variable that can be passed as an argument to more than one `branch` operator use the `branchCriteria` built-in method as shown below:

```
def criteria = branchCriteria {
    small: it < 10
    large: it > 10
}

Channel.from(1,2,30).branch(criteria).set { ch1 }
Channel.from(10,20,1).branch(criteria).set { ch2 }
```

### choice

### multiMap

### into

### separate

### tap
## Math operators

## Other operators