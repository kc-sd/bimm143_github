# Class 6: R Functions
Kyle Canturia (PID: A17502778)

- [Background](#background)
- [Our first function](#our-first-function)
- [A second function](#a-second-function)
- [A Protein generating function](#a-protein-generating-function)

## Background

All functions in R have at least 3 things:

- A **name** that we use to call the function.
- One or more input **arguments**
- The **body** the lines or R code that do the work

## Our first function

Let’s write a silly wee little function called `add()` to add some
numbers (the input arguments)

``` r
add <- function(x, y) {
  x + y 
}
```

Now we can use this function

``` r
add(100, 1)
```

    [1] 101

``` r
add(c(100, 1, 100), 1)
```

    [1] 101   2 101

> Q. What if I gave a multiple element vector to `x` and `y`?

``` r
add(x=c(100,1), y=c(100,1))
```

    [1] 200   2

> Q. What if I give three inputs to the function?

``` r
#add(x=c(100,1), y=1, z=1)
```

> Q. What if I give only one input to the add function?

``` r
addnew <- function(x, y=1) {
  x + y
}
```

``` r
addnew(x=100)
```

    [1] 101

``` r
addnew(c(100,1), 100)
```

    [1] 200 101

If we write our function with input arguments having no default value
then the user will be required to set them when they use the function.
We can give our input arguments”default” values by setting them equal to
some sensible value - e.g. y=1 in the `addnew()` function

## A second function

Let’s try something more interesting: Make a sequence generating tool..

The `sample()` function can be a useful starting point here:

``` r
sample(1:10, size=4)
```

    [1] 8 5 1 2

> Q. Generate 9 random numbers taken from the inpuut vector x=1:10

``` r
sample(1:10, size=9)
```

    [1]  3 10  1  4  8  9  6  2  7

> Q. Generate 12 random numbers taken from the inpuut vector x=1:10

``` r
sample(1:10, size=12, replace = T)
```

     [1] 3 2 6 8 8 4 3 4 6 9 3 1

> Q. Write code for the `sample()` function that generates nucleotides
> sequences of length 6

``` r
sample(x=c('a','t','g', 'c'), size=6, replace= T)
```

    [1] "t" "a" "g" "c" "g" "c"

> Q. Write a first function `generate_dna()` that returns a user
> specified length DNA sequence:

``` r
generate_dna <- function(length=6) {
  sample(c('A', 'T', 'G', 'C'), length, replace=T)
}
```

> **Key-Points** Every function in R looks fundamentally the same in
> terms of structure. Basically 3 things: name, input, and body

    name <- function(input) {
      body
    }

> Functions can have multiple inputs. These can be **required**
> arguments or **optional** arguments. With optional arguments having a
> set default value.

> Q. Modify and improve our `generate_dna()` function to return it’s
> generated sequence in a more standard format like “AGTAGTA” rather
> than the vector “A”,“C”,“G”,“A”

``` r
generate_dna <- function(length=6, fasta=T) {
  
  ans <- sample(c('A', 'T', 'G', 'C'), 
         length,
         replace=T)
  if(fasta) {
    cat("Single-element vector output")
    ans <- paste(ans, collapse = "")
  } else{
    cat("Multi-element vector output")
  }
  return(ans)
}

generate_dna(fasta=FALSE)
```

    Multi-element vector output

    [1] "T" "T" "T" "A" "C" "G"

``` r
generate_dna(fasta=T)
```

    Single-element vector output

    [1] "TTCCAG"

The `paste()` function - it’s job is to join up or stick together
(a.k.a. paste) input strings together

``` r
paste("alice", "loves R", sep=" ")
```

    [1] "alice loves R"

Flow control means where the R brain goes in your code

``` r
good_mood <- F

if(good_mood) {
  cat("Great!")
} else {
  cat("Bummer!")
}
```

    Bummer!

## A Protein generating function

> Q. Write a function that generates a user specified length protein
> sequence.

``` r
generate_protein <- function(length = 6, fasta=T) {
  if(length) {
     length = length
  }
  else{
    length <- sample(6:12, size=1,replace=T)
  }
  
  ans <- sample(c('G','A','P','V','L','I','M','F','Y','W','S','T','C','N','Q','K','R','H','D','E'), 
         length,
         replace=T)
  
  if(fasta) {
    ans <- paste(ans, collapse = "")
  } else{
    cat("")
  }
  return(ans)
}

generate_protein(fasta=T)
```

    [1] "DRIIIH"

> Q. Use that function to generate random protein sequences between
> length 6 and 12

``` r
generate_protein()
```

    [1] "EHIVFE"

``` r
for(i in 6:12) {
  # FASTA ID line ">id"
  cat(">", i, sep="", "\n")
  # Protein sequence line 
  cat(generate_protein(i), "\n")
}
```

    >6
    PWYLTW 
    >7
    ICNPNRR 
    >8
    GIRTDCNV 
    >9
    ENCHWMKEQ 
    >10
    YQFEFSTAVQ 
    >11
    QPEGAMGHIYN 
    >12
    HARFNQRMNGPC 

> Q. Are any of your sequences unique i.e. not found anywhere in nature?

Yes 9,10,11,12.
