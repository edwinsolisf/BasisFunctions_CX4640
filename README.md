---
Name: Edwin L. Solis F.
Topic: [26]
Title: Functional Analysis and the study of Basis Functions
---

Functional Analysis and the study of Basis Functions
=========================

Functional Analysis consists on the study of vector spaces and operations on that spaces that obey special kinds of properties defined to it. The most common example of such study subject is 
an n-dimensional real vector, that is an object that lives in $\mathbb{R}^n$ space. They obey certain rules related to dot products and norms while having operations done on them like linear maps
done through matrices. 

In a more general sense, it is possible to introduce an infinite dimensional complex vector space with inner product and norms, that is, a Hilbert space, in which vectors and functions become one and the same. With this vector behavior of functions, then all general theorems for vectors apply to functions providing advantages in their computation, representation, and operations. One of these advantages is the selection of special basis that can be used to represent any function by a linear combination of the basis just like in the case of vectors.
Such selection of basis functions should span the whole Hilbert Space.

As in linear algebra, there exist sets of orthogonal basis in which functions are orthogonal. For functions, orthogonality is defined in terms of the inner product being 0, in which these inner products are computed as integrals of the product of the two input functions.[^2] The usefulness of orthogonal polynomials comes from the ability to decompose any function into the linear combination of orthogonal functions and then utilize mathematical relations that can be apply in general to the orthogonal basis.

Orthogonal polynomials are used in many areas such as differential and integral equations, interpolation, numerical integration with Quadrature rules, eigenvalue problems, conformal mapping, and polynomial chaos in statistics.[^1][^3]

## Table of Contents
1. [History](#history)
2. [Mathematical Background](#mathematical-background)
3. [Basis Functions](#basis-functions)
4. [Orthogonal Functions](#orthogonal-functions)
    1. [Orthogonal Polynomials](#orthogonal-monomials)
    2. [Legendre Polynomials](#legendre-polynomials)
    3. [Leguerre Polynomials](#leguerre-polynomials)
    4. [Jacobi Polynomials](#jacobi-polynomials)
    5. [Other orthogonal functions](#other-orthogonal-functions)
        1. [Fourier Series](#fourier-series)
        2. [Rational Functions](#rational-functions)
5. [Numerical Analysis](#numerical-analysis)
6. [Examples](#examples)
7. [Further Reading](#further-reading)
8. [References](#references)

## History
While the idea of functions as quantities varying with respect to other parameters existed in the astronomical texts of ancient Greeks of the 100s BCE and in the works of the Persian mathematician Al-Biruni from the 11th century, the concrete and formal definition came from the development of Calculus from Leibniz, John Bernoulli, and Newton during the 19th Century.[^7]

Meanwhile, vectors appeared vaguely as in idea for grouping variables for solving systems of linear equations, mainly through the introduction of linear equations by Rene Descartes with geometry in 1637.[^8] It was not until 1844 that the core foundations of linear algebra took root from the publication "Theory of Extension" by Hermann Grassman. Here the concepts of linear combinations, basis, change of basis, subspaces, etc. are born. [^9]  Later this concepts would even be expanded to complex numbers such as the appearance of quaternions. However, the concept of vectors was completely based on linear algebra's matrices and linear equations, not an abstract malleable mathematical construct. Then, in 1888 Giuseppe Peano introduced the modern idea of abstract vector spaces. He proposed 4 concrete definitions that an _entity_ would follow, providing examples of such systems and even proposing the possibility of infinite dimensional vector spaces.[^4]

From Peano's realization, soon the idea of functions and infinite dimensional vector spaces would follow. In the 20th century the ideas of sequences and spaces of functions could be thought as vector spaces, yet the introduction of integral equations[^10], specifically two-square integrable functions that can have an inner product defined on an interval, showed its share of similarities with the Euclidean dot product. This ideas were presented by David Hilbert and Enrad Schmidt[^11]. Here the concept of functions as vector began to take shape as the definitions of orthogonality along with proofs of spectral decomposition. These key ideas would play a major role in defining a mathematical foundation in Quantum Mechanics, but its usefulness would travel far beyond from it.

By allowing functions and polynomials to be an element in a vector space, it allowed the possibilities for the selection of special sets of functions and the study of their properties. This way, the area of mathematics of Functional Analysis was born. Many of such special function sets have been compiled having different properties some of which we will present in this text.

## Mathematical Background
To understand in what sense functions are treatable as vectors, we present the definition of a _vector_, and proceed with some relations that may be useful for understanding the concept.

A vector space over a _number system_ $F$ (e.g. Integers $\mathbb{Z}$, Real $\mathbb{R}$, Complex $\mathbb{C}$, etc.) is a set $V$ that satisfies the following properties:[^12]

1. Associativity of Vector Addition: $\vec{u} + (\vec{v} + \vec{w}) = (\vec{u} + \vec{v}) + \vec{w}$
2. Communtativity of Vector Addition: $\vec{u} + \vec{v} = \vec{v} + \vec{u}$
3. Identity element of Vector Addition: $\vec{u} + \vec{0}= \vec{u}$
4. Inverse element of Vector Addition: $\vec{u} + (-\vec{u}) = \vec{0}$
5. Associativity of Scalar Multiplication: $(ab)\vec{u}=a(b\vec{u})$
6. Distributivity of Scalar Multiplication: $a(\vec{u} + \vec{v}) = a\vec{u} + a\vec{v}$ and $(a + b)\vec{u} = a\vec{u} + b\vec{u}$
7. Identity element of Scalar Multiplication: $1\vec{u} = \vec{u}$

where $\vec{u}$, $\vec{v}$ and $\vec{0}$ are elements in $V$ and $a,b$, and 1 are numbers in $F$.

These are the axioms, i.e. the definition or _statement assumed true without proof_, of vector spaces
Intuitively, for _geometrical_ vectors, this can be shown visually to be consistent to how vectors should behave.

Now, if we translate this concepts to functions, these properties still apply:
$(f+g)(x) = f(x) + g(x)$, $1f(x)=f(x)$, $a(f(x) + g(x)) = af(x) + ag(x)$, ... So, in essence, we can think of function and vectors as one and the same thing.

For the operation of __dot product__, more generally known as __inner product__, of a $n$-diemnsional vector with basis vectors $\hat{e_i}$ is defined as:
$$\vec{u}\cdot\vec{v}=\sum^n_{i=1}u_iv_i$$
with
$$u = \sum^n_{i=1}u_i\hat{e_i}$$
A function is a continous throughout is domain with infinitely many inputs, so a discrete definition with sumations will not work as it will diverge; therefore, the resolution occurs with defining the inner product with the continuous equivalent of the summation, the integral:

The inner product of two real functions $f$ and $g$ in the vector space $V$ over the interval $[a,b] is defined as:
$$
    \braket{f,g} = \int^b_af(x)\,g(x)\;\text{d}{x}
$$

More generally, functions can be complex and inner products can have specific weights depending on the topology:

$$
    \braket{f,g}_w = \int^b_af^*(x)\,g(x)\,w(x)\;\text{d}{x}
$$

where $f^*$ is the complex conjugate of $f$ and $w$ is the weigth of this particular inner product. Vector spaces that contain the inner product space are said to be __Inner Product Spaces__.

With this definition of inner products, it is also possible to define the norm of the function for the given interval. The __Euclidean Norm__, also called $L^2$ in the literature of Hilbert Spaces, for a function $f$ is defined as:
$$
    \big|\big|f\big|\big|_2 = \sqrt{\braket{f,f}_2} = \sqrt{\int^b_a\big|f(x)\big|^2\;\text{d}{x}}
$$

A function is said to be __normalized__ if its Euclidean Norm is equal to 1.

## Basis Functions

One of the main premises for the use of vectors is the ability to write a given vector as a linear combination of other vectors. These _other_ vectors are called __basis vectors__, and the vector space of the complete set of all possible combinations of those basis vectors is called the __span__.
If there is a set of $n$ basis vectors $e_i$ in vector space $V$, any vector $v$ in the subspace spanned by $e_i$, i.e. the span of $e_i$ can be written as:

$$\vec{v} = \sum^n_{i=1}v_i\vec{e_i}$$

The same principle applies to functions. A set function $S=\{\,f_i\;|\;1\leq i\leq n\}$ can represent any function in their span as a linear combination of those basis functions:

$$g(x) = \sum^n_{i=1}c_i\;f_i(x)$$

A simple example of this fact can be shown by a set of two basis functions: $f_0(x) = 1$ and $f_1(x)=x$. From the previous equation, it can be seen that any function $g(x)$ in the span of $\{f_0,f_1\}$ can be written as:
$$g(x) = c_0 + c_1 x$$
This is the equation of a line. Therefore, $g(x)$ can be any linear function or equivalently, the span of this set of basis function is the whole set of linear functions. Analyzing why this is the case, we can think of $c_0$ and $c_1$ as knobs. By only modifying $c_0$, $g(x)$ changes just in the intercept at $x=0$, but not the slope. Similarly, $c_1$ just modifies the slope, so the only linear function possible to construct with it are those that cross the origin. It is only when both can be modified is that the whole set of linear functions is attainable.

This analogy can be tried again with higher powers. If we add $f_2(x) = x^2$ to the set of basis functions, then $g(x)$ becomes:
$$g(x) = c_0 + c_1 x + c_2 x^2$$
The new knob $c_2$ opens the possibilities for $g(x)$ to be a second order polynomial, a quadratic.

By continuously adding more powers to this set, we construct what is called __Monomial Basis__. For a monomial basis with $n$ powers, the functions representable are given by the expression:
$$g(x) = \sum^{n-1}_{i=0}c_ix^i$$
which is just any polynomial of $n-1$ order.

In the case we let the number of monomials $n\to\infty$ for this set, then we would have:
$$g(x) = \sum^\infty_{i=0}c_ix^i$$
However, one must be careful to note that this expression must converge for it to be true. A set of examples of this would be Taylor expansions. For the exponential function, we would have:
$$g(x) = e^x = \sum^\infty_{n=0}\frac{1}{n!}x^n$$
where $c_n = \frac{1}{n!}$


## Orthogonal Functions



### Orthogonal Polynomials

### Legendre Polynomials 

### Leguerre Polynomials 

### Jacobi Polynomials

### Other orthogonal functions

#### Fourier Series[^6]

#### Rational Functions

#### Binary-valued functions

## Numerical Analysis

## Examples

## Further Reading

## References
[^1]: Gradimir V. Milovanovic, "Orthogonal Polynomial System and some Applications," Mathematical Institute of the Serbian Academy of Sciences and Arts, [Mathematical Institute SANU Website](http://www.mi.sanu.ac.rs/~gvm/radovi/inner.pdf), Accessed December 7, 2022.

[^2]: Gabor Szego, Gabor, "Orthogonal Polynomials," American Mathematical Society, Colloquium Publications, Volume 23, [Ohio State Universite Website](https://people.math.osu.edu/nevai.1/SZEGO/szego=szego1975=ops=OCR.pdf). Accessed December 7, 2022.

[^3]: T. J. Sullivan, "Orthogonal Polynomial Applications," Text in Applied Mathematics Book Series, Volume 63, [Springer Link](https://link.springer.com/chapter/10.1007/978-3-319-23395-6_8). Accessed December 7, 2022.

[^4]: Hubbert C. Kennedy, "Peano: Life and Works of Giuseppe Peano," pg. 23-25.  England, Reidel Publishing Company, 1980.

[^5]: The European Mathematical Society, "Hilbert Space," [Encylopedia of Mathematics](https://encyclopediaofmath.org/index.php?title=Hilbert_space). Accessed December 7, 2022

[^6]: Gerald B. Folland, "Fourier Analysis and Its Applications," The Sally Series, American Mathematical Society, 2009. 

[^7]: Medvedev F. A., "Scenes from the History of Real Functions," Birkhauser, 2012. Accessed from [Google Books](https://books.google.com/books?id=rqmABwAAQBAJ&dq=al+biruni+first+to+study+mathematical+function&pg=PA30#v=onepage&q=al%20biruni%20first%20to%20study%20mathematical%20function&f=false), December 7, 2022.

[^8]: Frank J. Swetz, "The Geometry of Rene Descartes," Mathematical Association of America. Accessed from [Mathemtatical Association of America Website](https://www.maa.org/publications/periodicals/convergence/the-geometry-of-rene-descartes), December 7, 2022.

[^9]: Desmond Fearnley-Sander, "Herman Grassman and the Creation of Linear Algebra," pg. 810-811 Mathematical Association of America, Accessed from [Mathematical Association of America Website](https://www.maa.org/sites/default/files/pdf/upload_library/22/Ford/DesmondFearnleySander.pdf), December 7, 2022.

[^10]: Nicholas Bourbaki, "Elements of Mathematics, General Topology", _Historical Note_, pg. 164-166, London, Springer, 1940

[^11]: Schmidt, E. Über die Auflösung linearer Gleichungen mit Unendlich vielen unbekannten. Rend. Circ. Matem. Palermo 25, 53–77 (1908). https://doi.org/10.1007/BF03029116

[^12]: Roman Steven, "Advanced Linear Algebra, Vector Spaces", Chapter 1, New York, Springer, 2005.