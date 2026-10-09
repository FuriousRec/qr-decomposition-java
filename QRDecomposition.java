public class QRDecomposition {
    private static final double EPSILON = 1e-10;

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
     * Performs the QR decomposition using the Gram-Schmidt process.
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

            // Remove the parts of the vector that point in previous Q directions.
            for (int previousColumn = 0; previousColumn < currentColumn; previousColumn++) {

                // Get a previous column of Q.
                double[] previousQColumn = getColumn(Q, previousColumn);

                // Calculate how much the current input column points in that Q direction.
                R[previousColumn][currentColumn] =
                        dotProduct(previousQColumn, getColumn(inputMatrix, currentColumn));

                // Subtract that part from the remaining vector.
                for (int row = 0; row < rowCount; row++) {
                    remainingVector[row] =
                            remainingVector[row] - R[previousColumn][currentColumn] * previousQColumn[row];
                }
            }

            // Store the length of the remaining vector in R.
            R[currentColumn][currentColumn] = norm(remainingVector);

            // If length is almost 0, the columns are linearly dependent.
            if (Math.abs(R[currentColumn][currentColumn]) < EPSILON) {
                throw new IllegalArgumentException("Matrix columns are linearly dependent.");
            }

            // Normalize the remaining vector and store it as the current column of Q.
            for (int row = 0; row < rowCount; row++) {
                Q[row][currentColumn] = remainingVector[row] / R[currentColumn][currentColumn];
            }
        }

        return new QRResult(Q, R);
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

    // Calculates length of a vector.
    private static double norm(double[] vector) {
        return Math.sqrt(dotProduct(vector, vector));
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
