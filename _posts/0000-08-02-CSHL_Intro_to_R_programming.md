---
title: "Introduction to R"
date: '0000-08-02'
feature_image: assets/genvis-dna-bg_optimized_v1a.png
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
categories: Module-00-Setup
---
**Author : Katie Campbell, UCLA**<br>

***
# Introduction to R

This session will cover the basics of R programming, from setting up your environment to basic data analysis. In this tutorial, you will:

-   Become comfortable navigating your filesystem and environment from the R console
-   Understanding, reading, and manipulating data structures
-   Change or summarize datasets for basic analysis

------------------------------------------------------------------------

## Introduction to R programming

R is a powerful programming language and software environment used for statistical computing, data analysis, and graphical representation of data. It is widely used among statisticians, data analysts, and researchers for data mining and statistical software development.

### Prerequisite: Files for today's session

Today's session will utilize two input files, which we will use for practice. These can be downloaded into your instance with the following commands:

``` bash
curl https://raw.githubusercontent.com/ParkerICI/MORRISON-1-public/refs/heads/main/RNASeq/RNA-CancerCell-MORRISON1-metadata.tsv > intro_r_metadata.tsv
curl https://raw.githubusercontent.com/ParkerICI/MORRISON-1-public/refs/heads/main/RNASeq/data/RNA-CancerCell-MORRISON1-combat_batch_corrected-logcpm-all_samples.tsv.zip > intro_r_dataset.tsv.zip
```

The downloaded `intro_r_metadata.tsv` file contains annotation of a set of RNAseq samples from patients with melanoma treated with immunotherapies. The `intro_r_dataset.tsv.zip` file contains batch effect-corrected gene expression values for all of the samples in this dataset.

### Why use R?

-   **Open Source**: R is free to use and open-source.
-   **Extensive Packages**: Thousands of packages available for various statistical techniques.
-   **Strong Community Support**: Active community contributes to continuous improvement.
-   **Cross-Platform**: Works on Windows, macOS, and Linux.

### Getting started

We will launch R in our instance and be programming within the terminal directly. This would be equivalent to writing code in the "Console" panel of R studio. Launch R using the command:

``` bash
R
```

The terminal should print out the following information:

``` bash
R version 4.4.2 (2024-10-31) -- "Pile of Leaves"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 
```

The *\>* character indicates that you are now within the R environment.

After the course, you can download R from [CRAN](https://cran.r-project.org/) on your personal computer. Once R is installed, you can run it from your terminal or using an Integrated Development Environment (IDE) like [RStudio](https://www.rstudio.com/), which makes it convenient to code, debug, and analyze data within a convenient user interface.

Some of the commands in this tutorial have a comment character (`#`) in front of them. These are commands that you don't have to run in today's session, but can come back to on your own time. Commenting your code in scripts will also support the documentation and clarity of your scripts.

------------------------------------------------------------------------

## Working directory and packages

### Working directory

The working directory is the folder where R reads and saves files by default.

You can check your working directory by:

``` r
getwd()
```

If you are interested in moving to a different directory (for example, where your files are being stored), you can change the working directory:

``` r
# setwd("/new/path/")
```

You can replace `"/new/path"` with your desired path. Note that directory paths are subtly different between Mac and Windows.

### Packages

Packages extend R's functionality by providing additional functions and datasets. First, they need to be installed before they can be loaded. Many of the packages used in this course will already be available in your instance.

``` r
# install.packages("tidyverse")
```

Once packages are installed, they can be loaded using *library(package_name)*. This is usually the first thing that you run in your script or console to load all of the functions you will be using. Tidyverse is a collection of R packages designed for data science, including packages like `ggplot2` for data visualization and `dplyr` for data manipulation.

``` r
library(tidyverse)
```

As you learn commands, you can include *?* before the command in order to see a description of what the command does, the parameters and inputs, and the structure of the output. For example, you can see the documentation for *list.files()* with the following command:

``` r
?list.files()
# Type q to exit the help documentation
```

This shows that *list.files()* will list all of the files in a directory or folder as a character vector. You can change the directory to search in with the parameter *path* or specify a type of file with *pattern*. As you perform data analysis in R over the course of this workshop, this is a helpful way to explore what commands are doing and how you can use them.

*Try this:* What files are in your current working directory?

You can also check the files in the parent directory, subdirectories, or based upon file path patterns:

``` r
list.files(path = "..")
list.files(pattern = ".tsv")
```

------------------------------------------------------------------------

## Variables and the environment

In python, `=` assigns the value on the right to the name of the variable on the left In R, `<-` or `=` can be used to assign a value on the right to the name of the variable on the left

``` r
age <- 32
first_name <- 'Katie'
```

### Displaying values

`print()` can be used to show a value, but you can also just type the variable name.

``` r
print(age)
age
print(first_name)
first_name
```

The value 32 is shown without quotes, while "Katie" is shown with quotes, since 32 is an **integer** (or **double**) and "Katie" is a **character** string. One other **variable type** is a **logical** (TRUE or FALSE value). Variable types can be determined using the commands `typeof()`.

``` r
typeof(age)
typeof(32)

typeof(first_name)

typeof(TRUE)
```

Note that if we want to print out multiple variables, we can't simply list them in `print()` the way that we did in python.

``` r
print(first_name, "is", age, "years old")
```

Instead, we can use `paste()` function to combine this statement:

``` r
print(paste(first_name, "is", age, "years old"))
```

`paste()` will concatenate the values, separating them by spaces. *Try this:* What happens if you use the `paste0()` function instead?

As in python, variables must be created before they are used.

``` r
print(last_name)
```

This results in the error: *Error: object 'last_name' not found*

### The global environment

The environment is where R stores variables and functions you've created. In RStudio, all of your stored variables are visibly listed in the upper right panel. You can also list the objects of your environment with the command `ls()`:

``` r
ls()
```

As we store more variables, this list grows:

``` r
last_name <- "Campbell"
ls()
```

We can also remove variables from our environment:

``` r
rm(age)
ls()
```

Or we can remove all objects and clear the environment entirely:

``` r
rm(list = ls())
```

------------------------------------------------------------------------

## Creating data structures

Data structures can be combined into one- and two-dimensional objects in R, including vectors, lists, matrices, and data frames.

-   **Vectors** are created using the `c()` function, and all of the elements must be of the same type.
-   **Lists** are created using the `list()` function and can contain elements of multiple types.
-   **Matrices** are created using the `matrix()` function and are two-dimensional arrays that contain one type.
-   **Data frames** are created using the `data.frame()` function, are two-dimensional, and columns contain vectors of the same type.

### Vectors and lists

You can create a numeric vector of only numbers. We can use the function `str()` to further explore the dimensions and types in the object.

``` r
numbers <- c(1, 2, 3, 4, 5)
typeof(numbers)
str(numbers)
```

A numeric vector can also be generated by a sequence of numbers:

``` r
sequence <- 1:10
sequence

skipping_sequence <- seq(from = 1, to = 10, by = 2)
skipping_sequence

repeated_sequence <- rep(x = 1, times = 10)
repeated_sequence
```

A character vector only contains character strings:

``` r
names <- c("Alice", "Bob", "Charlie", "Daniel", "Edwin")
typeof(names)
str(names)

names_sequence <- rep(x = "Your Name", 10)
names_sequence
```

Note that the structure command `str()` shows not only the type of vector, but also the length of it.

If you try to create a vector of different types, R will automatically coerce the variable into the most general type. For example, if you try to create a vector of numbers and letters, it will coerce the entire vector into characters:

``` r
values <- c(1, 2, 3, 4, "Alice", "Bob")
values
str(values)

more_values <- c(numbers, names)
str(more_values)
```

Notice that all of the values print out with quotation marks, and `str()` shows that it is a character (chr) vector.

### Lists

Lists are unique from vectors

### Matrices

Matrices are generated using the `matrix()` function. You can create a matrix by feeding the command a specific set of values. Look at the parameters for `matrix()`:

``` r
?matrix()
```

The default of `matrix()` is to create a matrix with one column. To organize these values into a specific size, you need to fill in the parameters `nrow` and `ncol`:

``` r
matrix(13:24) # Creates a 12x1 matrix of these values

# These are all the same
matrix(13:24, nrow = 3)
matrix(13:24, ncol = 4)
matrix(13:24, nrow = 3, ncol = 4)

# These fill in the values by row, rather than in the order of the columns
matrix(13:24, nrow = 3, byrow = TRUE)
matrix(13:24, ncol = 4, byrow = TRUE)
matrix(13:24, nrow = 3, ncol = 4, byrow = TRUE)
```

### Data frames

Data frames are used for storing two-dimensional data structures, where each row represents a series of observations. Each column represents a variable, which is a vector of one type.

``` r
df <- data.frame(
  name = names,
  age = numbers,
  flower = c(rep("rose", 3), rep("petunia", 2)),
  favorite_color = c('red','blue','green','yellow','macaroni and cheese'),
  has_dog = c(TRUE, FALSE, TRUE, FALSE, FALSE)
)
```

You can view your data frame using `head()`:

``` r
head(df)
head(df, 3)
```

Vectors can be extracted from the data frame using `data_frame$column_name`:

``` r
df$name
df$age
```

Some useful commands to understand your data structures, in general, include summarizing the rows, columns, and names of these dimensions:

``` r
str(df)
nrow(df) # Number of rows
ncol(df) # Number of columns
dim(df) # Dimensions of the data frame (row, column)
colnames(df)
names(df)
```

Data frames store the type of each variable as a vector:

``` r
summary(df)
```

#### Adding, removing, and modifying columns

``` r
df$new_column <- "my new column"
df

df$new_column <- NULL
df

df$lower_name <- tolower(df$name)
df

df$has_dog <- ifelse(df$has_dog, "yes", "no")
df

df$has_dog <- ifelse(df$has_dog == "yes", TRUE, FALSE)
df
```

In `tidyverse`, you can use the `mutate()` column to do the same thing. You just have to overwrite the variable each time:

```r
mutate(df, new_column = "my new column")
df

df <- mutate(df, new_column = "my new column")
df

# You can also include a list of these changes in one
mutate(df, new_column = "change it again", new_name = toupper(name))
```

------------------------------------------------------------------------

## Indexing

Unlike python, R uses one-based indexing. So the index of the first element is 1, not 0.

### Indexing one dimensional objects: Vectors, lists

For one-dimensional objects, like vectors or lists, we just use brackets to extract the specific index:

``` r
numbers[1]

df$name[3]
```

This is also how we overwrite these values:

``` r
df$name[3] <- "Frederick"
df
```

``` r
colnames(df) # outputs a vector of the column names
colnames(df)[2] <- "age_in_years"
colnames(df)
```

Lists are sometimes named or can be more complex, since not every element of the list has to be the same structure or type. This is where double brackets may come into play:

``` r
complicated_list <- list("first" = c(1,2,3,4), "second" = df, "third" = 1087.29)
str(complicated_list)
complicated_list
complicated_list[[1]]
complicated_list[['first']]
complicated_list$first
complicated_list[['first']][3]
```

### Indexes of two dimensional objects: Arrays, data frames

For two dimensional objects, you have to specify both the rows and columns that you're interested in viewing, using `[row, column]`:

``` r
# This can be done with indices
my_matrix <- matrix(13:24, nrow = 3)
my_matrix
my_matrix[1,4]

# Or if the rows/columns are named
colnames(my_matrix) <- LETTERS[1:4]
row.names(my_matrix) <- paste0("row", 1:3)
my_matrix
my_matrix['row1','C']

# To get all of the rows or all of the columns you leave the field blank
# BUT you still need the comma
my_matrix['row3',]

# Without the comma...
my_matrix['A']
my_matrix['row2']
```

------------------------------------------------------------------------

## String Manipulation

The `stringr` package in R (one of the packages in `tidyverse`) simplifies these tasks with easy-to-use functions that can handle typical string operations.

### Finding patterns

Finding specific sequences or motifs within biological sequences is a common task.

``` r
sequence <- "ATGCGTACGTTGACA"
motif <- "CGT"
str_locate(sequence, motif)
```

### Replacing substrings

Modifying sequences by replacing specific nucleotides or amino acids.

``` r
dna_sequence <- "ATGCGTACGTTGACT"
rna_sequence <- str_replace_all(dna_sequence, "T", "U")
print(rna_sequence)
```

### Substring extraction

Extracting parts of sequences, such as cutting out genes or regions of interest.

``` r
extracted_sequence <- str_sub(sequence, 3, 8)
print(extracted_sequence)
```

### Length calculation

Determining the length of sequences.

``` r
sequence_length <- str_length(sequence)
print(sequence_length)
```

### Case conversion

Converting uppercase to lowercase, or vice versa.

``` r
sequence_upper <- str_to_upper(sequence)
print(sequence_upper)
```

## Splitting strings

Splitting sequences into arrays, useful for reading fasta files or analyzing codons.

``` r
codons <- str_sub(sequence, seq(1, str_length(sequence), by = 3), seq(3, str_length(sequence), by = 3))
print(codons)
```

*Try this:* What if our sequence length wasn't a multiple of three?

## Counting specific characters

Counting occurrences of specific nucleotides or amino acids.

``` r
guanine_count <- str_count(sequence, "G")
print(guanine_count)
```

------------------------------------------------------------------------

## Logical operations

Data analysis will require filtering data types using comparison operators:

-   `==` checks to see whether two values are equivalent. A common error is to only use one `=`, which is used to set variables in functions.
-   `!=` checks if two values are not equivalent
-   `>`, `<`, `>=`, and `<=` check for greater/less than (or equal to) for comparing numbers

You can use these operators to subset vectors:

``` r
numbers<3
which(numbers<3)
numbers[numbers<3]
numbers[which(numbers<3)]
```

Conditional statements can be combined using logical operators:

-   `&` requires that both statements are true
-   `|` requires that one of the statements is true
-   `!` requires that the opposite of the statement is true

``` r
numbers>1 & numbers<5
which(numbers>1 | numbers<5)
numbers[numbers>1 & numbers<5]
numbers[which(numbers>1 & numbers<5)]
```

When we want to use these operations to filter data frames, we can specify the rows from the data frame by filtering the vectors stored in the columns. Note that we need to include the comma to indicate that we want all of the rows of this subsetted data frame:

``` r
df[df$Age < 30, ]
```

Alternatively, we can also use the `filter()` command from the `dplyr` package:

``` r
filter(df, Age < 30)
```

------------------------------------------------------------------------

## Reading in your own data

Instead of typing out all of the data elements, you can import data from various file formats. Some file types are easily read in by `base R`. Note that there are packages that may be helpful that are specific to certain input file types.

### Reading tab- or comma-delimited files

It is important to know the format of the file you are trying to read. The extension *tsv* indicates that the values are tab-separated.

``` r
metadata <- read.table("intro_r_metadata.tsv")
```

What is wrong with this data when you evaluate it?

``` r
head(metadata)
```

You'll notice that the column names are in the first row of data, and all of the columns are labeled as "V" and the number corresponding to the column. It helps to be specific about how R should read in the data:

``` r
metadata <- read.table("intro_r_metadata.tsv", sep = "\t", head = TRUE)
head(metadata)
str(metadata)
```

### Reading Excel Files

The `readxl` package can be used to read in Excel files:

``` r
# install.packages("readxl")
# library(readxl)
```

``` r
# excel_data <- read_excel("path/to/excel.xlsx", sheet = "Sheet1")
```

### Exploring the data

You can get a glimpse of the data by only looking at the first few lines:

``` r
head(metadata)
```

You can quickly perform summary statistics on a dataset or a vector:

``` r
summary(metadata)
summary(metadata$subject.age)
```

**Try this:**

-   Filter data to the samples that are cutaneous, and save this as a new data frame.
-   What is the average age of the patients in this filtered dataset?

### Factors

Currently, all of the character strings have been read in as characters. This is helpful, but oftentimes our data is more structured than that. For example, the `bor` column in metadata encodes clinical outcomes, and our brains think of these as an ordinal variable, not just a categorical one. That is, we think of the outcomes in order from best to worst: CR > PR > SD > PD (Complete response is better than partial response, which is better than stable disease, then progressive disease). 

A `factor` encodes this information using `levels`. When we read in data, strings are read in as character vectors, but we can coerce them to become factors for our data summary:

```r
summary(metadata$bor)
unique(metadata$bor)

metadata$bor <- factor(metadata$bor)
summary(metadata$bor)

metadata$bor <- factor(metadata$bor, levels = c("CR","PR","SD","PD"))
summary(metadata$bor)
```

**Try this:**

-   Filter metadata to the samples from patients that had complete response to therapy, and save this as a new data frame.
-   How many of these tumors came from cutaneous tumors?


### Sorting data

Like numbers, factors have a logical order to them. So when we enforce the data to have this particular order, we can also arrange the data in that order.

```r
numbers
order(numbers)
order(-numbers)

arrange(metadata, bor)
arrange(metadata, age)
arrange(metadata, -age)


arrange(metadata, -bor) # GIVES A WARNING
arrange(metadata, desc(bor)) # Better for factors
```

------------------------------------------------------------------------

## Complex data manipulation

### Long and Wide Data Formats

Long and wide data formats are two common ways of structuring data, each with its own advantages and use cases.

### Long Format

In the long format, also known as the "tidy" format, each observation is represented by a single row in the dataset. This format is characterized by having:

-   Multiple rows, each corresponding to a single observation or measurement.
-   One column for the variable being measured.
-   Additional columns to store metadata or grouping variables.

**Advantages**:

-   Facilitates easy analysis and manipulation, especially when using tools like Tidyverse packages in R.
-   Suitable for data that follow the "one observation per row" principle, such as time series or longitudinal data.

### Wide Format

In the wide format, each observation is represented by a single row, but with multiple columns corresponding to different variables. This format is characterized by:

-   One row per observation.
-   Each variable is represented by a separate column.

**Advantages**:

-   Can be easier to understand for simple datasets with few variables.
-   May be more convenient for certain types of analysis or visualization.

### Choosing Between Long and Wide Formats

The choice between long and wide formats depends on factors such as the nature of the data, the analysis tasks, and personal preference. Long format is often preferred for its flexibility and compatibility with modern data analysis tools, while wide format may be suitable for simpler datasets or specific analysis requirements.

#### Long to Wide

``` r
library(tidyr)

# Example long format data
long_data <- data.frame(
  Subject = c("A", "A", "B", "B"),
  Time = c(1, 2, 1, 2),
  Measurement = c(10, 15, 12, 18)
)

# Convert long format data to wide format
wide_data <- spread(long_data, key = Time, value = Measurement)

# View the wide format data
print(wide_data)
```

#### Wide to Long

``` r
library(tidyr)

# Example wide format data
wide_data <- data.frame(
  Subject = c("A", "B"),
  Time1 = c(10, 12),
  Time2 = c(15, 18)
)

# Convert wide format data to long format
long_data <- gather(wide_data, key = Time, value = Measurement, -Subject)

# View the long format data
print(long_data)
```

#### Example: Gene expression data

Let's work with a real dataset: `intro_r_dataset.tsv.zip`. This file is compressed, so we can't read it in by the default `read.table` function. This is a circumstance where I use the `data.table` package and the `fread` function, which automatically recognizes the compression format to read in the file:

```r
read.table("intro_r_dataset.tsv.zip") # DOESN'T WORK

# install.packages('data.table')
data <- fread("intro_r_dataset.tsv.zip")

dim(data)
colnames(data)
str(data)
```

*Try this:* Create a long data frame, where the key is called "sample.id" and the value column contains the gene expression ("cpm"). Call this new dataset "data_long".

### Merging Data

Merging allows combining data from different sources. This is common in analyzing biological data. Joins and merging are common operations used to combine multiple datasets based on common variables or keys. In Tidyverse, these operations are typically performed using functions from the `dplyr` package.

#### Types of Joins:

##### Inner Join (`inner_join()`):

An inner join combines rows from two datasets where there is a match based on a common key, retaining only the rows with matching keys from both datasets.

##### Left Join (`left_join()`):

A left join combines all rows from the first (left) dataset with matching rows from the second (right) dataset based on a common key. If there is no match in the second dataset, missing values are filled in.

##### Right Join (`right_join()`):

Similar to a left join, but it retains all rows from the second (right) dataset and fills in missing values for non-matching rows from the first (left) dataset.

##### Full Join (`full_join()`):

A full join combines all rows from both datasets, filling in missing values where there are no matches.

##### Semi-Join (`semi_join()`):

A semi-join returns only rows from the first dataset where there are matching rows in the second dataset, based on a common key.

##### Anti-Join (`anti_join()`):

An anti-join returns only rows from the first dataset that do not have matching rows in the second dataset, based on a common key.

##### Merge (`merge()`):

The `merge()` function is a base R function used to merge datasets based on common columns or keys. It performs similar operations to joins in `dplyr`, but with slightly different syntax and behavior.

##### Example:

``` r
library(dplyr)

# Example datasets
df1 <- data.frame(ID = c(1, 2, 3), Name = c("Alice", "Bob", "Charlie"))
df2 <- data.frame(ID = c(2, 3, 4), Score = c(85, 90, 95))

# Inner join
inner_merged <- inner_join(df1, df2, by = "ID")

# Left join
left_merged <- left_join(df1, df2, by = "ID")

# Right join
right_merged <- right_join(df1, df2, by = "ID")

# Full join
full_merged <- full_join(df1, df2, by = "ID")

# Semi-join
semi_merged <- semi_join(df1, df2, by = "ID")

# Anti-join
anti_merged <- anti_join(df1, df2, by = "ID")
```

#### Example: Merging our metadata with our gene expression data

Now let's merge our `data_long` with `metadata` so that our gene expression data also contains our sample annotation.

```r
left_join(data_long, metadata)

```

------------------------------------------------------------------------

## Chaining commands, groupby(), and summarise()

So far, we've used individual commands to accomplish several tasks, but sometimes we want to do multiple things in one line of code. The `%>%` is called a pipe and is used to chain commands together.

```r
df
df %>% 
  mutate(age_in_days = age_in_years*365) %>% 
  filter(age < 500)
  
metadata %>%
  filter(sample.tumor.type == "cutaneous") %>%
  arrange(bor)
```

`groupby()` allows us to apply individual functions to grouped objects, and we can use `summarise()` to perform functions within those individual groups:

```r
metadata %>%
  group_by(sample.tumor.type, bor) %>% 
  summarise(total = n(), mean_age = mean(subject.age, na.rm = TRUE)) # n() performs a count
  
metadata %>%
  group_by(sample.tumor.type, bor) %>% count() # This is another way to quickly count groups
```

*Advanced exercise:* Create a wide data frame that summarizes the total number of responders (defined by `bor`) per tumor type (in each row).

------------------------------------------------------------------------

## Repeating tasks/functions

The `groupby()` functionality in tidyverse allows you to perform many metrics across individual groups, so you don't have to create a filtered dataset over and over again. Base R also has a series of functions to make it easier to repeat the same function over and over again.

### `apply`

The apply() function in R is a powerful tool for applying a function to the rows or columns of a matrix or data frame. It is particularly useful for performing operations across a dataset without needing to write explicit loops. The syntax for apply() is:


```r
apply(X, margin, function, ...)

# X: This is the array or matrix on which you want to apply the function.
# margin: A value that specifies whether to apply the function over rows (1), columns (2), or both (c(1, 2)).
# function: The function you want to apply to each row or column.
```

To calculate the sum of each row in a matrix:
```r
# Create a matrix
my_matrix <- matrix(1:9, nrow=3)

# Apply sum function across rows
row_sums <- apply(my_matrix, 1, sum)
print(row_sums)
```

To find the mean of each column in a data frame:
```r
# Create a data frame
df <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))

# Apply mean function across columns
column_means <- apply(df, 2, mean)
print(column_means)
```

### `sapply` and `lappy`

- `lapply()` returns a list, regardless of the output of each application of the function.
- `sapply()` attempts to simplify the result into a vector or matrix if possible. If simplification is not possible, it returns a list similar to `lapply()`.

Suppose you have a list of numerical vectors and you want to compute the sum of each vector. Here's how you could use lapply():
```R
# Define a list of vectors
num_list <- list(c(1, 2, 3), c(4, 5), c(6, 7, 8, 9))

# Use lapply to apply the sum function
list_sums <- lapply(num_list, sum)
print(list_sums)
```

Using the same list of numerical vectors, if you use sapply() to compute the sum, the function will try to simplify the output into a vector:

```r
# Use sapply to apply the sum function
vector_sums <- sapply(num_list, sum)
print(vector_sums)
```

When to Use Each
 
 - `lapply()`: When you need the robustness of a list output, especially when dealing with heterogeneous data or when the function can return variable lengths or types.
 - `sapply()`: When you are working with homogeneous data and prefer a simplified output such as a vector or matrix, assuming the lengths and types are consistent across elements.

### Advanced: Using enframe() and unnest() to create a data frame

Another option, if you enjoy using tidyverse and groupby() is to convert a list of objects into a data structure amenable to this:

```r
data1 <- data.frame("subject" = paste0("subject", sample(1:10, 5, replace = FALSE)),
  "age" = sample(1:100, 5))
data2 <- data.frame("subject" = paste0("subject", sample(11:20, 5, replace = FALSE)),
  "age" = sample(1:100, 5))
data3 <- data.frame("subject" = paste0("subject", sample(21:30, 5, replace = FALSE)),
  "age" = sample(1:100, 5))

list_of_data <- list("data1" = data1, "data2" = data2, "data3" = data3)

enframe(list_of_data)

enframe(list_of_data) %>% unnest(value)
```

`enframe()` will convert a list of objects into a tibble (a type of data frame), and each list element is encoded in the value column.

`unnest(value)` will unnest the data frames from the `value` column into individual rows.

This is very helpful if you ever have many files that are the same type, and you're trying to create one master data frame.

```r
# my_files <- list.files(path = "/path/to/data", pattern = "expression.tsv")
# merged_data <- lapply(my_files, fread) %>% enframe() %>% unnest(value)
```

------------------------------------------------------------------------

## Additional exercises

Try these exercises, using the `metadata` and `data_long` objects that you read into R:

* What sample had the highest expression of *B2M*?
* Which `cohort` had the highest average expression of *B2M*?
* Which gene had the highest average expression across all patients?
* Which response group (`bor`) had the highest average expression of *NLRC5*?
* What is the median expression of *HLA-A* within each tumor type/response group?

------------------------------------------------------------------------

## Save R objects for future use

Throughout the course, we will complete one stage of analysis and save objects from the environment to individual files for future use. If you want to export all objects in the environment, you can use `save.image("path/to/file.rds")`. Alternatively, individual files are stored and then loaded later using the following commands:

``` r
saveRDS(data, file = "testdata.rds")

# Next time you open R, you can reload the object with:
load("testdata.rds")
```

------------------------------------------------------------------------

## Closing R

When you want to quit R in your terminal, you can type in the following commmand:

``` r
q()
```

The console will ask you:

``` r
Save workspace image? [y/n/c]: 
```

Answering "yes" and saving your workspace image will create a hidden file called *.Rdata* in your current working directory. This will store all of the data objects in your global environment to that .Rdata file and can be automatically loaded for future use. If you open R again in the future, this will automatically reload all of the stored data objects, recreating the environment. This is helpful for continuing analysis from where you previously left off.

Answering "no" will not save your current environment, and you will need to rerun all of the code that you previously ran. Oftentimes, all of your code will be stored in an R script and will reproduce your prior analysis.

------------------------------------------------------------------------

## Additional Resources

-   [R Documentation](https://www.rdocumentation.org/)
-   [R for Data Science](https://r4ds.had.co.nz/) by Garrett Grolemund and Hadley Wickham
-   [Advanced R](https://adv-r.hadley.nz/) by Hadley Wickham


