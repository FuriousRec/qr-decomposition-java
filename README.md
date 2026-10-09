# QR Decomposition

This project implements QR decomposition in Java using the Gram-Schmidt process.

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
column. Columns whose residual norm is below the fixed `1e-10` threshold are
also rejected.

The tests cover reconstruction, all three non-finite values in every position
of a 3-by-2 matrix, invalid shapes, dependent columns, and rejection of
non-finite test results. GitHub Actions compiles and runs them on Java 8 and 17.

This is an educational classical Gram-Schmidt implementation. Finite input
validation does not guarantee numerical stability: very large finite values
can still overflow, and the absolute dependency threshold is scale-sensitive.
