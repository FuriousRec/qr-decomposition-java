public class QRDecomposition {
    private static final double RANK_TOLERANCE = 1e-12;

    /**
     * Helper class to store the two resulting matrices.
     * This is static so it can be used without creating an instance of the outer class.
     */
    static class QRResult {
        double[][] Q;
        double[][] R;

        QRResult(double[][] Q, double[][] R) {
            this.Q = Q;
            this.R = R;
        }
    }

    /**
     * Performs the QR decomposition using the modified Gram-Schmidt process.
     */
    public static QRResult qrDecomposition(double[][] inputMatrix) {
        validateInput(inputMatrix);

        int rowCount = inputMatrix.length;
        int columnCount = inputMatrix[0].length;

        // Creating new matrices for us to store after we decompose our main matrix.
        double[][] Q = new double[rowCount][columnCount];  // Q stores the orthonormal vectors.
        double[][] R = new double[columnCount][columnCount];  // R stores values used to rebuild the original matrix.

        // Go through every column of the input matrix.
        for (int currentColumn = 0; currentColumn < columnCount; currentColumn++) {
            double[] remainingVector = getColumn(inputMatrix, currentColumn);
            double originalColumnNorm = norm(remainingVector);

            if (!Double.isFinite(originalColumnNorm)) {
                throw new IllegalArgumentException(
                        "Matrix column norm exceeds the finite range (column "
                                + currentColumn + ").");
            }

            // Remove the parts of the vector that point in previous Q directions.
            for (int previousColumn = 0; previousColumn < currentColumn; previousColumn++) {

                // Get a previous column of Q.
                double[] previousQColumn = getColumn(Q, previousColumn);

                // Project the already-reduced vector. This is modified Gram-Schmidt:
                // using the residual improves orthogonality for nearly dependent columns.
                R[previousColumn][currentColumn] =
                        dotProduct(previousQColumn, remainingVector);

                // Subtract that part from the remaining vector.
                for (int row = 0; row < rowCount; row++) {
                    remainingVector[row] =
                            remainingVector[row] - R[previousColumn][currentColumn] * previousQColumn[row];
                }
            }

            // Store the length of the remaining vector in R.
            double residualNorm = norm(remainingVector);
            R[currentColumn][currentColumn] = residualNorm;

            // Compare relative to this column's original scale. An absolute cutoff
            // incorrectly rejects small but independent matrices.
            if (residualNorm <= RANK_TOLERANCE * originalColumnNorm) {
                throw new IllegalArgumentException("Matrix columns are linearly dependent.");
            }

            // Normalize the remaining vector and store it as the current column of Q.
            for (int row = 0; row < rowCount; row++) {
                Q[row][currentColumn] = remainingVector[row] / residualNorm;
            }
        }

        return new QRResult(Q, R);
    }

    /**
     * Solves the full-column-rank least-squares problem min ||Ax - b||_2.
     *
     * <p>For square matrices this is an ordinary linear solve. For tall
     * matrices it returns the vector whose residual is orthogonal to every
     * column of A.</p>
     */
    public static double[] solveLeastSquares(
            double[][] inputMatrix,
            double[] rightHandSide) {

        validateInput(inputMatrix);
        validateRightHandSide(rightHandSide, inputMatrix.length);

        QRResult decomposition = qrDecomposition(inputMatrix);
        int columnCount = inputMatrix[0].length;
        double[] projected = new double[columnCount];

        // y = Q^T b
        for (int column = 0; column < columnCount; column++) {
            projected[column] =
                    dotProduct(getColumn(decomposition.Q, column), rightHandSide);
            requireFiniteLeastSquaresValue(projected[column]);
        }

        // Solve the upper-triangular system Rx = y by back substitution.
        double[] solution = new double[columnCount];
        for (int row = columnCount - 1; row >= 0; row--) {
            double remaining = projected[row];

            for (int column = row + 1; column < columnCount; column++) {
                remaining -= decomposition.R[row][column] * solution[column];
            }

            solution[row] = remaining / decomposition.R[row][row];
            requireFiniteLeastSquaresValue(solution[row]);
        }

        return solution;
    }

    /**
     * Solves several full-column-rank least-squares problems with one QR
     * decomposition.
     *
     * <p>Each column of {@code rightHandSides} is an independent right-hand
     * side. If A is m-by-n and B is m-by-k, the returned solution is n-by-k
     * and minimizes ||AX-B|| for every column of B.</p>
     */
    public static double[][] solveLeastSquares(
            double[][] inputMatrix,
            double[][] rightHandSides) {

        validateInput(inputMatrix);
        validateRightHandSides(rightHandSides, inputMatrix.length);

        QRResult decomposition = qrDecomposition(inputMatrix);
        int rowCount = inputMatrix.length;
        int columnCount = inputMatrix[0].length;
        int rightHandSideCount = rightHandSides[0].length;
        double[][] projected =
                new double[columnCount][rightHandSideCount];

        // Y = Q^T B. The QR decomposition is shared by every right-hand side.
        for (int column = 0; column < columnCount; column++) {
            for (int rightHandSide = 0;
                    rightHandSide < rightHandSideCount;
                    rightHandSide++) {

                double value = 0.0;
                for (int row = 0; row < rowCount; row++) {
                    value += decomposition.Q[row][column]
                            * rightHandSides[row][rightHandSide];
                }
                requireFiniteLeastSquaresValue(value);
                projected[column][rightHandSide] = value;
            }
        }

        double[][] solution =
                new double[columnCount][rightHandSideCount];

        // Solve RX = Y by back substitution, one result column at a time.
        for (int rightHandSide = 0;
                rightHandSide < rightHandSideCount;
                rightHandSide++) {

            for (int row = columnCount - 1; row >= 0; row--) {
                double remaining = projected[row][rightHandSide];

                for (int column = row + 1;
                        column < columnCount;
                        column++) {

                    remaining -= decomposition.R[row][column]
                            * solution[column][rightHandSide];
                }

                solution[row][rightHandSide] =
                        remaining / decomposition.R[row][row];
                requireFiniteLeastSquaresValue(
                        solution[row][rightHandSide]);
            }
        }

        return solution;
    }

    private static void validateRightHandSides(
            double[][] rightHandSides,
            int expectedRows) {

        if (rightHandSides == null
                || rightHandSides.length == 0
                || rightHandSides[0] == null
                || rightHandSides[0].length == 0) {

            throw new IllegalArgumentException(
                    "Right-hand side matrix must not be empty.");
        }
        if (rightHandSides.length != expectedRows) {
            throw new IllegalArgumentException(
                    "Right-hand side matrix row count must equal "
                            + "the matrix row count.");
        }

        int columns = rightHandSides[0].length;
        for (int row = 0; row < rightHandSides.length; row++) {
            if (rightHandSides[row] == null
                    || rightHandSides[row].length != columns) {

                throw new IllegalArgumentException(
                        "Right-hand side matrix must be rectangular.");
            }
            for (int column = 0; column < columns; column++) {
                if (!Double.isFinite(rightHandSides[row][column])) {
                    throw new IllegalArgumentException(
                            "Right-hand side matrix entries must be finite "
                                    + "(row " + row
                                    + ", column " + column + ").");
                }
            }
        }
    }

    private static void validateRightHandSide(
            double[] rightHandSide,
            int expectedLength) {

        if (rightHandSide == null) {
            throw new IllegalArgumentException(
                    "Right-hand side must not be null.");
        }
        if (rightHandSide.length != expectedLength) {
            throw new IllegalArgumentException(
                    "Right-hand side length must equal the matrix row count.");
        }
        for (int index = 0; index < rightHandSide.length; index++) {
            if (!Double.isFinite(rightHandSide[index])) {
                throw new IllegalArgumentException(
                        "Right-hand side entries must be finite (index "
                                + index + ").");
            }
        }
    }

    private static void requireFiniteLeastSquaresValue(double value) {
        if (!Double.isFinite(value)) {
            throw new ArithmeticException(
                    "Least-squares computation produced a non-finite result.");
        }
    }

    private static void validateInput(double[][] matrix) {
        if (matrix == null || matrix.length == 0 || matrix[0] == null || matrix[0].length == 0) {
            throw new IllegalArgumentException("Matrix must not be empty.");
        }

        int columns = matrix[0].length;

        if (matrix.length < columns) {
            throw new IllegalArgumentException("QR decomposition requires at least as many rows as columns.");
        }

        for (int row = 0; row < matrix.length; row++) {
            if (matrix[row] == null || matrix[row].length != columns) {
                throw new IllegalArgumentException("Matrix must be rectangular.");
            }
            for (int column = 0; column < columns; column++) {
                if (!Double.isFinite(matrix[row][column])) {
                    throw new IllegalArgumentException(
                            "Matrix entries must be finite (row " + row + ", column " + column + ").");
                }
            }
        }
    }

    // --- Helper Methods ---
    private static double[] getColumn(double[][] matrix, int columnIndex) {
        double[] column = new double[matrix.length];

        for (int row = 0; row < matrix.length; row++) {
            column[row] = matrix[row][columnIndex];
        }

        return column;
    }

    private static double dotProduct(double[] firstVector, double[] secondVector) {
        double result = 0.0;

        for (int index = 0; index < firstVector.length; index++) {
            result += firstVector[index] * secondVector[index];
        }

        return result;
    }

    // Calculates length without squaring first, avoiding needless overflow/underflow.
    private static double norm(double[] vector) {
        double result = 0.0;
        for (double value : vector) {
            result = Math.hypot(result, value);
        }
        return result;
    }

    // Prints a matrix nicely.
    private static void printMatrix(String name, double[][] matrix) {
        System.out.println(name + ":");

        for (int row = 0; row < matrix.length; row++) {
            for (int column = 0; column < matrix[row].length; column++) {
                System.out.printf("%10.5f ", matrix[row][column]);
            }
            System.out.println();
        }

        System.out.println();
    }

    public static void main(String[] args) {
        double[][] exampleMatrix = {
                {1, 1},
                {1, 0},
                {0, 1}
        };

        QRResult result = qrDecomposition(exampleMatrix);

        printMatrix("Q", result.Q);
        printMatrix("R", result.R);
    }
}

