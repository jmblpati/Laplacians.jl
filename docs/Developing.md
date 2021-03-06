# Developing Laplacians.jl

[TOC]

## Learn to use git

* If you don't know anything about git, then just know that you should make a branch for you own code.  Type

~~~
git checkout -b MyName
~~~

* Now, read about Git.  I recommend the book
[Pro Git](https://git-scm.com/book/en/v2), which is available online for free.

* Stop thinking about Git like subversion or dropbox.

* The master branch will be the one for public consumption. It should (mostly) work.

* You should also read the 
[section of the Julia docs about building packages.](http://docs.julialang.org/en/release-0.4/manual/packages/#package-development)

## Fast code?

Just go for it.
Don't worry about writing fast code at first.
Just get it to work.
We can speed it up later.


Within some of the files, I am keeping old, unoptimized versions of code around for comparison (and for satisfaction).  I will give them the name "XSlow"

## Documentation

This documentation is still very rough.
It is generated by a combination of Markdown and semi-automatic generation.  The steps to generate and improve it are:

* Edit Markdown files in the `docs` directory.  For example, you could use MacDown to do this.
* If you want to add a new page to the documention, create one.  Edit the file mkdocs.yml so show where it should appear.
* Add docstrings to everything that needs it, and in particular to the routines you create.  The API is built from the docstrings.  To build the API, type

~~~julia
include("docs/build.jl")
~~~

or run `julia docs/build.jl` from the root directory.

* Run `mkdocs build` in the root directory to regenerate the documentation from the Markdown.

* Once you like the documentation, you can upload it with 

~~~
mkdocs gh-deploy --clean -b gh-pages
~~~

* *Warning:* mkdocs deletes everying in gh-pages that it does not put there itself.

* If you create a Julia notebook that you would like to include as documentation.   You should
   put it in the notebooks directory (.julia/v0.4/Laplacians/notebooks) and then link to it's page on GitHub.  While it seems that one should convert it to html (and one can), and then include it in MkDocs, MkDocs does something funny to the resulting html that does not look nice.




## Parametric Types

A sparse matrix has two types associated with it: the types of its indices (some sort of integer) and the types of its values (some sort of number).  Most of the code has been written so that once these types are fixed, the type of everything else in the function has been too.  This is accomplished by putting curly braces after a function name, with the names of the types that we want to use in the braces.  For example,

~~~julia
shortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, start::Ti)
~~~

`Tv`, sometimes written `Tval` denotes the types of the values, and `Ti` or `Tind` denotes the types of the indices.  This function will only be called if the node from which we compute the shortest paths, `start` is of type `Ti`.  Inside the code, whenever we write something like `pArray = zeros(Ti,n)`, it creates an array of zeros of type Ti.  Using these parameteric types is *much* faster than leaving the types unfixed.

### Data structures:

* `IntHeap` a heap that stores small integers (like indices of nodes in a graph) and that makes deletion fast.  Was much faster than using Julia's more general heap.

### Interface issue:
There are many different sorts of things that our code could be passing around.  For example, kruskal returns a graph as a sparse matrix.  But, we could use a format that is more specialized for trees, like the RootedTree type.  At some point, when we optimize code, we will need to figure out the right interfaces between routines.  For example, some routines symmetrize at the end.  This is slow, and should be skipped if not necessary.  It also doubles storage.

### Writing tests:
I haven't written any yet.  I'll admit that I'm using the notebooks as tests.  If I can run all the cells, then it's all good.

## Integration with other packages.

There are other graph packages that we might want to sometimes use.

* [Graphs.jl](http://github.com/JuliaLang/Graphs.jl) : I found this one to be too slow and awkward to be useful.
* [LightGraphs.jl](http://github.com/JuliaGraphs/LightGraphs.jl) : this looks more promising.  We will have to check it out.  It is reasonably fast, and the code looks pretty.

