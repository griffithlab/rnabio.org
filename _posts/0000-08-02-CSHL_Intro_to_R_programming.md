---
output: pdf_document
---
# Introduction to R

This session will cover the basics of R programming, from setting up your environment to basic data analysis. In this tutorial, you will:

-   Become comfortable navigating your filesystem and environment from the R console
-   Understandin, reading, and manipulating data structures
-   Change or summarize datasets for basic analysis and plotting

------------------------------------------------------------------------

## Introduction to R Programming {#introduction-to-r-programming}

R is a powerful programming language and software environment used for statistical computing, data analysis, and graphical representation of data. It is widely used among statisticians, data analysts, and researchers for data mining and statistical software development.

### Why Use R?

-   **Open Source**: R is free to use and open-source.
-   **Extensive Packages**: Thousands of packages available for various statistical techniques.
-   **Strong Community Support**: Active community contributes to continuous improvement.
-   **Cross-Platform**: Works on Windows, macOS, and Linux.

### Getting started

R can be launched from bash, simply by running the command:

``` bash
R
```

This will automatically open the R console, printing the following message:

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

The *\>* character indicates that you are now within the R environment. As you learn commands, you can include *?* before the command in order to see a description of what the command does, the parameters and inputs, and the structure of the output. For example, you can see the documentation for *list.files()* with the following command:

``` r
?list.files()
```

This shows that *list.files()* will list all of the files in a directory or folder as a character vector. You can change the directory to search in with the paramter *path* or specify a type of file with *pattern*. As you perform data analysis in R over the course of this workshop, this is a helpful way to explore what commands are doing and how you can use them.

After the course, you can download R from [CRAN](https://cran.r-project.org/) on your personal computer. Once R is installed, you can run it from your terminal or using an Integrated Development Environment (IDE) like [RStudio](https://www.rstudio.com/), which makes it convenient to code, debug, and analyze data within a convenient user interface.

------------------------------------------------------------------------

## Working Directory, Packages, and the Global Environment {#working-directory-packages-and-the-global-environment}

### Working Directory

The working directory is the folder where R reads and saves files by default.

You can check your working directory by:

``` r
getwd()
```

If you are interested in moving to a different directory (for example, where your files are being stored), you can change the working directory:

``` r
setwd("/home/ubuntu/workspace")
```

*You can replace `"/home/ubuntu/workspace"` with your desired path.*

### Packages

Packages extend R's functionality by providing additional functions and datasets. First, they need to be installed before they can be loaded. Many of the packages used in this course will already be available in your instance.

``` r
install.packages("tidyverse")
```

Once packages are installed, they can be loaded using *library(package_name)*. This is usually the first thing that you run in your script or console to load all of the functions you will be using. Tidyverse is a collection of R packages designed for data science, including packages like `ggplot2` for data visualization and `dplyr` for data manipulation.

``` r
library(tidyverse)
```

### Global Environment

The global environment is where R stores variables and functions you've created. To see all of the objects in your global environment, use:

``` r
ls()
```

As you create objects, you can also remove them:

``` r
new_number <- 1
ls()
rm(new_number)
```

To remove all objects, clearing your global environment, you can run:

``` r
rm(list = ls())
```

------------------------------------------------------------------------

## Variables {#variables}

Variables store data values. In R, you can assign values to variables using `<-` or `=` to store them in the global environment.

``` r
x <- 10
y = 5
```

The most common types of variables are:

-   **Numeric**: Integers and doubles (e.g., `x <- 10`)
-   **Character**: Text strings (e.g., `name <- "John"`)
-   **Logical**: TRUE or FALSE values (e.g., `is_valid <- TRUE`)

You can view variables by printing them to the console:

``` r
print(x)
```

Or simply type the variable name:

``` r
x
```

------------------------------------------------------------------------

## Creating Data Structures {#creating-data}

You can create various data structures in R, such as vectors, matrices, lists, and data frames.

### Vectors

A vector is a sequence of data elements of the same basic type.

A numeric vector only contains numbers:

``` r
numbers <- c(1, 2, 3, 4, 5)
```

A character vector only contains character strings:

``` r
names <- c("Alice", "Bob", "Charlie")
```

If you try to list of different types, R will default to the most general type. For example, if you try to create a vector of numbers and letters, it will coerce the entire vector into characters:

``` r
values <- c(1, 2, 3, 4, "Alice", "Bob")
```

The commands *str()* and *typeof()* can help you understand the size and type of data you are analyzing:

``` r
str(numbers)
str(names)
str(values)

typeof(numbers)
typeof(names)
typeof(values)
```

Specific values can be extracted from a vector by specifying the index:

``` r
length(numbers)
numbers[3]
```

**Try this:**

-   What is the 4th value in the vector `values`?
-   How long is the vector `names`?
-   What are the first 3 values in `names`?

### Number sequences

A vector containing a sequence of numbers can be generated by:

``` r
sequence <- 1:10
```

Or using the `seq()` function:

``` r
sequence <- seq(from = 1, to = 10, by = 2)
```

### Repetition

Repeat elements using `rep()`.

``` r
repeated <- rep(x = 1, times = 5)
```

### Matrices

Matrices are two-dimensional arrays that contain

``` r
matrix_data <- matrix(1:9, nrow = 3, ncol = 3)
```

### Lists

Unlike vectors, lists can contain elements of different types.

``` r
my_list <- list(numbers = numbers, names = names)
```

### Data Frames

Data frames are used for storing two-dimensional data structures, where each row represents a series of observations. Each column represents a variable

``` r
df <- data.frame(
  ID = 1:3,
  Name = c("Alice", "Bob", "Charlie"),
  Age = c(25, 30, 35)
)
```

The variables from a data frame can be accessed by:

``` r
df$ID
```

To further explore the variables and structure of the data set, you can use the following commands:

``` r
colnames(df)
str(df)
nrow(df)
ncol(df)
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

## Reading In and Interpreting Data {#reading-in-and-interpreting-data}

Instead of typing out all of the data elements, you can import data from various file formats like CSV, Excel, and others.

### Reading tab- or comma-delimited files

It is important to know the format of the file you are trying to read. The extension *tsv* indicates that the values are tab-separated.

``` r
data <- read.table("/home/ubuntu/workspace/rnaseq/de/deseq2/DE_all_genes_DESeq2.tsv")
```

What is wrong with this data when you evaluate it?

``` r
head(data)
```

You'll notice that the column names are in the first row of data, and all of the columns are labeled as "V" and the number corresponding to the column. It helps to be specific about how R should read in the data:

``` r
data <- read.table("/home/ubuntu/workspace/rnaseq/de/deseq2/DE_all_genes_DESeq2.tsv", sep = "\t", head = TRUE)
head(data)
str(data)
```

### Reading Excel Files

The `readxl` package can be used to read in Excel files:

``` r
install.packages("readxl")
library(readxl)
```

``` r
excel_data <- read_excel("path/to/excel.xlsx", sheet = "Sheet1")
```

### Viewing Data

You can get a glimpse of the data by only looking at the first few lines:

``` r
head(data)
```

You can quickly perform summary statistics on a dataset or a vector:

``` r
summary(data)
summary(data$HBR_Rep2)
summary(numbers)
```

**Try this:**

-   Filter data to the genes with an adjusted p value cutoff of 0.05
-   How many genes are significantly different in this dataset?
-   How many genes are significantly upregulated in this dataset?

------------------------------------------------------------------------

## Manipulating Data Frames {#manipulating-data-frames}

Data manipulation is essential for preparing data for analysis.

### Selecting Columns

``` r
# Select by column names
data_subset <- data[c("Symbol", "log2FoldChange","pvalue")]

# Select by column index
data_subset <- data[c(7,3,5)]

# Using dplyr
data_subset <- select(data, Symbol, log2FoldChange, pvalue)
```

### Filtering Rows

``` r
filtered_data <- filter(data, log2FoldChange>2)
```

### Adding New Columns

``` r
# Add a new column directly
data$my_new_column <- "this is my new column"

# Use mutate() from dplyr
# data <- mutate(data, my_new_column = "this is my new column")
```

### Modifying Existing Columns

``` r
# Convert all pvalues to -log10pvalue
data$pvalue <- -log10(data$pvalue)

# mutate() also modifies existing columns
# data <- mutate(data, pvalue = -log10(pvalue))
```

### Removing Columns

``` r
# Remove 'lfcSE' column
data$lfcSE <- NULL
```

### Sorting Data

``` r
# Sort by pvalue in ascending order (remember this is converted to the -log10 pvalue)
sorted_data <- data[order(data$pvalue), ]

# Sort by pvalue in descending order
sorted_data <- data[order(-data$pvalue), ]

# dplyr has the `arrange()` function for sorting values
arrange(data, -pvalue)
```

**Try this:**

-   Create a new data object called `de_result` by reading in the file `/home/ubuntu/workspace/rnaseq/de/deseq2/DE_sig_genes_DESeq2.tsv`
-   Filter the data to the genes that are significant with an adjusted p value cutoff of 0.01 and a log2-fold change of 2
-   What is the median log2FoldChange of these genes? How many genes are there?
-   When you sort the genes by gene name, what is the first one that comes up alphabetically?

------------------------------------------------------------------------

## Basic Plotting and Statistics {#basic-plotting-and-statistics}

## Long and Wide Data Formats

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

## Long to Wide

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

## Wide to Long

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

# Merging Data

Merging allows combining data from different sources. This is common in analyzing biological data.

## Joins and Merging of Data in Tidyverse

Joins and merging are common operations used to combine multiple datasets based on common variables or keys. In Tidyverse, these operations are typically performed using functions from the `dplyr` package.

### Types of Joins:

#### Inner Join (`inner_join()`):

An inner join combines rows from two datasets where there is a match based on a common key, retaining only the rows with matching keys from both datasets.

#### Left Join (`left_join()`):

A left join combines all rows from the first (left) dataset with matching rows from the second (right) dataset based on a common key. If there is no match in the second dataset, missing values are filled in.

#### Right Join (`right_join()`):

Similar to a left join, but it retains all rows from the second (right) dataset and fills in missing values for non-matching rows from the first (left) dataset.

#### Full Join (`full_join()`):

A full join combines all rows from both datasets, filling in missing values where there are no matches.

#### Semi-Join (`semi_join()`):

A semi-join returns only rows from the first dataset where there are matching rows in the second dataset, based on a common key.

#### Anti-Join (`anti_join()`):

An anti-join returns only rows from the first dataset that do not have matching rows in the second dataset, based on a common key.

### Merging Data:

#### Merge (`merge()`):

The `merge()` function is a base R function used to merge datasets based on common columns or keys. It performs similar operations to joins in `dplyr`, but with slightly different syntax and behavior.

### Example:

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

# Plotting in ggplot2

The core idea behind `ggplot2` is the concept of a "grammar of graphics". This concept provides a systematic way to describe and build graphical presentations such as charts and plots. The grammar itself is a set of independent components that can be composed in many different ways. This grammar includes elements like:

-   Data: The raw data that you want to visualize.
-   Aesthetics (`aes`): Defines how data are mapped to color, size, shape, and other visual properties.
-   Geometries (`geom`): The geometric objects in a plot—lines, points, bars, etc.
-   Scales: Transformations applied to data before it is visualized, including scales for colors, sizes, and shapes.
-   Coordinate systems: The space in which data is plotted.
-   Facets: Used for creating plots with multiple panels (small multiple plots).
-   Statistical transformations (stat): Summary statistics that can be applied to data before it is visualized, such as counting or averaging.
-   Themes: Visual styles and layout configurations for the plot.

Here’s how you generally use ggplot2 to create a plot:

-   Start with `ggplot()`: Set up the data and, optionally, define default mappings between variables and their aesthetics.
-   Add layers: Add layers to the plot using geom\_ functions, such as `geom_point()` for scatter plots, `geom_line()` for line graphs, and so on.\
-   Adjust the scales: Customize the scales used for aesthetics such as color, size, and x-y coordinates.
-   Modify the coordinate system: Choose a coordinate system.
-   Add facets: If necessary, add facets to create a multi-panel plot.
-   Apply a theme: Customize the appearance of the plot using themes.

``` r
library(ggplot2)

# Sample data
df <- data.frame(
  x = rnorm(100),
  y = rnorm(100),
  group = factor(rep(1:2, each = 50))
)

# Creating a scatter plot
p <- ggplot(df, aes(x = x, y = y, color = group)) + 
  geom_point() +
  theme_minimal() +
  labs(title = "Scatter Plot Example", x = "X Axis", y = "Y Axis")

# Saving the scatter plot
pdf("my_file.pdf")
p
dev.off()
```

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
