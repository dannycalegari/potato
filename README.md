potato
======

Program for "potato packing", after Bobenko-Pinkall-Springborn. This is an ingredient in the Hurwitz realization
problem (details elsewhere . . .)

The program accepts files in several different formats. The default has the following format:

vertices

valence(0) A(0,0), L(0,0) A(0,1) L(0,1) . . A(0,valence(0)-1), L(0,valence(0)-1)

valence(1) A(1,0), . . .

. . .

valence(vertices-1) A(vertices-1,0), . . .


where vertices is an integer and specifies the number of vertices, valence(i) is an integer and specifies the valence of
vertex i, A(i,j) is an integer and specifies which vertex is the jth neighbor of vertex i (in cyclic order), and L(i,j) is
a dbl such that e^{L/2} is the length of the edge from i to A(i,j). Here "dbl" is defined in the preprocessor, and is by
default a long double.

A variation on the default formats, called "no length format", omits the Ls and sets them all to 0 (so that all edge 
lengths have unadjusted lengths 1).

A third variation is "dual format", produced by Laurent Bartholdi's program. In this format a file is in the form

VERTICES n

infty

FACES m

i j k l_jk l_ik l_ij

... (repeated m times)


Meaning: there are n vertices; vertex infty (in {1...n}) is to be put at infinity; there are m faces, all triangular; each
line describes a triangle. The line above describes the triangle with vertices i j k on its CCW boundary, and lengths of edges
jk, ik, ij (*not* L as above!). The words VERTICES and FACES must appear as text in the file in the appropriate places.

To run the program on a file in default format, run ./potato filename

To run the program on a file in no length format, run ./potato -nl filename

To run the program on a file in dual format, run ./potato -df filename

For a file in dual format, the program outputs a list of packed centers stereographically projected to the unit sphere in 
R^3, and saves this to a file "center_list.txt". It also draws an eps output of the Euclidean packing (minus the star of the
"infinite" vertex) to a file "potato_packing.eps". An excerpt from sample output is:
![sample output](https://raw.github.com/dannycalegari/potato/master/sample_potato_packing.pdf)
