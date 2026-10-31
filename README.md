# QR Decomposition

This project implements QR decomposition in Java using the modified Gram-Schmidt process.

QR decomposition splits a matrix `A` into two matrices:

```text
A = Q * R
```

where `Q` contains orthonormal columns and `R` is upper triangular.

## Files

```text
QRDecomposition.java
QRDecompositionTest.java
```

## Least-squares solving

The same implementation can solve square systems and full-column-rank
overdetermined systems without forming the less stable normal equations:

```java
double[][] design = {
    {1, 0},
    {1, 1},
    {1, 2}
};
double[] observations = {1, 2, 2};

double[] coefficients =
        QRDecomposition.solveLeastSquares(design, observations);
// coefficients: intercept = 7/6, slope = 1/2
```

The solver calculates `Q^T b` and then uses back substitution on `R`.
For overdetermined inputs it minimizes the Euclidean residual `||Ax-b||`.

Multiple output vectors can be solved together without recomputing the QR
decomposition. Each column of `observations` is an independent right-hand
side, and each result column contains its coefficients:

```java
double[][] observations = {
    {1, 2},
    {4, 1},
    {7, 0}
};

double[][] coefficients =
        QRDecomposition.solveLeastSquares(design, observations);
// columns: [1, 3] and [2, -1]
```

## Run the Demo

```bash
javac --release 8 QRDecomposition.java
java QRDecomposition
```

## Run the Tests

```bash
javac --release 8 QRDecomposition.java QRDecompositionTest.java
java QRDecompositionTest
```

Expected output:

```text
All QR regression tests passed.
```

## Input validation

Input must be a nonempty rectangular matrix with at least as many rows as
columns. Every entry must be finite: `NaN` and positive/negative infinity are
rejected with an `IllegalArgumentException` identifying the zero-based row and
column. Rank detection compares each residual norm with its original column
norm using a relative `1e-12` tolerance, so uniformly scaling a matrix does not
change whether its columns are considered dependent.

The tests cover reconstruction, orthogonality on an ill-conditioned matrix,
all three non-finite values in every position of a 3-by-2 matrix, invalid
shapes, dependent columns, extreme finite scales, vector and multi-output least-squares solutions,
normal-equation residuals, invalid right-hand sides, and rejection of non-finite
test results. GitHub Actions compiles and runs them on Java 8 and 17.

This is an educational modified Gram-Schmidt implementation. Each projection
uses the progressively reduced residual vector, which retains orthogonality
better than classical Gram-Schmidt for nearly dependent columns. It is still
not as robust as Householder QR. The norm calculation uses `Math.hypot` to avoid
overflow and underflow from squaring finite values; a column whose true norm
exceeds the finite `double` range is rejected with a clear error.

