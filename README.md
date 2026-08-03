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
Test passed.
```
