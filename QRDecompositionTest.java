public class QRDecompositionTest {
    private static final double TOLERANCE = 1e-9;

    public static void main(String[] args) {
        testQrRebuildsOriginalMatrix();
        testIllConditionedMatrixKeepsOrthonormalColumns();
        testScaleExtremesRemainStable();
        testSolvesExactSquareSystem();
        testFitsOverdeterminedData();
        testLeastSquaresResidualIsOrthogonal();
        testLeastSquaresScaleExtremes();
        testInvalidRightHandSidesAreRejected();
        testNonFiniteEntriesAreRejected();
        testInvalidShapesAreRejected();
        testDependentColumnsAreRejected();
        testFiniteAssertionRejectsNaN();
        System.out.println("All QR regression tests passed.");
    }

    private static void testQrRebuildsOriginalMatrix() {
        double[][] originalMatrix = {
                {1, 1},
                {1, 0},
                {0, 1}
        };

        QRDecomposition.QRResult result = QRDecomposition.qrDecomposition(originalMatrix);
        double[][] rebuiltMatrix = multiply(result.Q, result.R);

        assertMatrixClose(originalMatrix, rebuiltMatrix);
    }

    private static void testIllConditionedMatrixKeepsOrthonormalColumns() {
        double epsilon = 1e-6;
        double[][] illConditionedMatrix = {
                {1, 1, 1},
                {epsilon, 0, 0},
                {0, epsilon, 0},
                {0, 0, epsilon}
        };

        QRDecomposition.QRResult result =
                QRDecomposition.qrDecomposition(illConditionedMatrix);

        // A classical Gram-Schmidt projection loses about 5e-5 orthogonality
        // on this matrix. Modified Gram-Schmidt keeps the error below 1e-8.
        for (int firstColumn = 0; firstColumn < result.Q[0].length; firstColumn++) {
            for (int secondColumn = 0;
                    secondColumn < result.Q[0].length;
                    secondColumn++) {

                double expected = firstColumn == secondColumn ? 1.0 : 0.0;
                double actual = columnDotProduct(
                        result.Q, firstColumn, secondColumn);
                if (!Double.isFinite(actual)
                        || Math.abs(expected - actual) > 1e-8) {
                    throw new AssertionError(
                            "Q columns are not orthonormal: expected "
                                    + expected + " but got " + actual);
                }
            }
        }

        assertMatrixClose(
                illConditionedMatrix,
                multiply(result.Q, result.R));
    }

    private static void testScaleExtremesRemainStable() {
        for (double scale : new double[] {1e-200, 1e200}) {
            double[][] matrix = {
                    {scale, scale},
                    {scale, -scale},
                    {0, scale}
            };

            QRDecomposition.QRResult result =
                    QRDecomposition.qrDecomposition(matrix);

            assertOrthonormal(result.Q, 1e-12);
            assertMatrixRelativeClose(matrix, multiply(result.Q, result.R), 1e-12);
        }
    }

    private static void testSolvesExactSquareSystem() {
        double[][] matrix = {
                {2, 1},
                {1, 3}
        };
        double[] rightHandSide = {1, 2};

        assertVectorClose(
                new double[] {0.2, 0.6},
                QRDecomposition.solveLeastSquares(matrix, rightHandSide),
                1e-12);
    }

    private static void testFitsOverdeterminedData() {
        // Best-fit line through (0,1), (1,2), (2,2):
        // intercept = 7/6 and slope = 1/2.
        double[][] designMatrix = {
                {1, 0},
                {1, 1},
                {1, 2}
        };
        double[] observations = {1, 2, 2};

        assertVectorClose(
                new double[] {7.0 / 6.0, 0.5},
                QRDecomposition.solveLeastSquares(
                        designMatrix, observations),
                1e-12);
    }

    private static void testLeastSquaresResidualIsOrthogonal() {
        double[][] matrix = {
                {1, 0},
                {1, 1},
                {1, 2},
                {1, 3}
        };
        double[] rightHandSide = {1, 2, 2, 5};
        double[] solution =
                QRDecomposition.solveLeastSquares(matrix, rightHandSide);

        for (int column = 0; column < matrix[0].length; column++) {
            double normalEquationResidual = 0.0;
            for (int row = 0; row < matrix.length; row++) {
                double fitted = 0.0;
                for (int inner = 0; inner < solution.length; inner++) {
                    fitted += matrix[row][inner] * solution[inner];
                }
                normalEquationResidual +=
                        matrix[row][column]
                                * (fitted - rightHandSide[row]);
            }

            if (Math.abs(normalEquationResidual) > 1e-12) {
                throw new AssertionError(
                        "Residual is not orthogonal to column "
                                + column + ": " + normalEquationResidual);
            }
        }
    }

    private static void testLeastSquaresScaleExtremes() {
        for (double scale : new double[] {1e-200, 1e200}) {
            double[][] matrix = {
                    {scale, 0},
                    {0, scale},
                    {scale, scale}
            };
            double[] rightHandSide = {
                    2 * scale,
                    -3 * scale,
                    -scale
            };

            assertVectorClose(
                    new double[] {2, -3},
                    QRDecomposition.solveLeastSquares(
                            matrix, rightHandSide),
                    1e-12);
        }
    }

    private static void testInvalidRightHandSidesAreRejected() {
        double[][] matrix = {
                {1, 0},
                {0, 1},
                {1, 1}
        };

        assertLeastSquaresRejected(
                matrix,
                null,
                "Right-hand side must not be null.");
        assertLeastSquaresRejected(
                matrix,
                new double[] {1, 2},
                "Right-hand side length must equal the matrix row count.");
        assertLeastSquaresRejected(
                matrix,
                new double[] {1, Double.NaN, 3},
                "Right-hand side entries must be finite (index 1).");
        assertLeastSquaresRejected(
                matrix,
                new double[] {1, Double.POSITIVE_INFINITY, 3},
                "Right-hand side entries must be finite (index 1).");
    }

    private static void assertLeastSquaresRejected(
            double[][] matrix,
            double[] rightHandSide,
            String expectedMessage) {

        try {
            QRDecomposition.solveLeastSquares(matrix, rightHandSide);
        } catch (IllegalArgumentException expected) {
            if (!expectedMessage.equals(expected.getMessage())) {
                throw new AssertionError(
                        "Unexpected validation message: "
                                + expected.getMessage());
            }
            return;
        }
        throw new AssertionError(
                "Expected invalid right-hand side to be rejected.");
    }

    private static void testNonFiniteEntriesAreRejected() {
        double[] invalidValues = {Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
        for (double invalid : invalidValues) {
            // Exercise every position, including entries beyond the first row/column.
            for (int row = 0; row < 3; row++) {
                for (int column = 0; column < 2; column++) {
                    double[][] matrix = {{1, 1}, {1, 0}, {0, 1}};
                    matrix[row][column] = invalid;
                    assertRejected(matrix, "Matrix entries must be finite (row "
                            + row + ", column " + column + ").");
                }
            }
        }
    }

    private static void testInvalidShapesAreRejected() {
        assertRejected(null, "Matrix must not be empty.");
        assertRejected(new double[0][], "Matrix must not be empty.");
        assertRejected(new double[][] {null}, "Matrix must not be empty.");
        assertRejected(new double[][] {{}}, "Matrix must not be empty.");
        assertRejected(new double[][] {{1, 2}}, 
                "QR decomposition requires at least as many rows as columns.");
        assertRejected(new double[][] {{1}, null}, "Matrix must be rectangular.");
        assertRejected(new double[][] {{1}, {2, 3}}, "Matrix must be rectangular.");
    }

    private static void testDependentColumnsAreRejected() {
        assertRejected(new double[][] {{1, 2}, {2, 4}},
                "Matrix columns are linearly dependent.");
    }

    private static void assertRejected(double[][] matrix, String expectedMessage) {
        try {
            QRDecomposition.qrDecomposition(matrix);
        } catch (IllegalArgumentException expected) {
            if (!expectedMessage.equals(expected.getMessage())) {
                throw new AssertionError("Unexpected validation message: " + expected.getMessage());
            }
            return;
        }
        throw new AssertionError("Expected invalid matrix to be rejected.");
    }

    private static void testFiniteAssertionRejectsNaN() {
        for (double invalid : new double[] {Double.NaN, Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY}) {
            try {
                assertMatrixClose(new double[][] {{1}}, new double[][] {{invalid}});
            } catch (AssertionError expected) {
                continue;
            }
            throw new AssertionError("Comparison accepted a non-finite result.");
        }
    }

    private static double columnDotProduct(
            double[][] matrix, int firstColumn, int secondColumn) {

        double result = 0.0;
        for (double[] row : matrix) {
            result += row[firstColumn] * row[secondColumn];
        }
        return result;
    }

    private static void assertOrthonormal(double[][] matrix, double tolerance) {
        for (int firstColumn = 0; firstColumn < matrix[0].length; firstColumn++) {
            for (int secondColumn = 0;
                    secondColumn < matrix[0].length;
                    secondColumn++) {

                double expected = firstColumn == secondColumn ? 1.0 : 0.0;
                double actual = columnDotProduct(matrix, firstColumn, secondColumn);
                if (!Double.isFinite(actual)
                        || Math.abs(expected - actual) > tolerance) {
                    throw new AssertionError(
                            "Q columns are not orthonormal: expected "
                                    + expected + " but got " + actual);
                }
            }
        }
    }

    private static double[][] multiply(double[][] left, double[][] right) {
        double[][] result = new double[left.length][right[0].length];

        for (int row = 0; row < left.length; row++) {
            for (int column = 0; column < right[0].length; column++) {
                for (int sharedIndex = 0; sharedIndex < right.length; sharedIndex++) {
                    result[row][column] += left[row][sharedIndex] * right[sharedIndex][column];
                }
            }
        }

        return result;
    }

    private static void assertVectorClose(
            double[] expected,
            double[] actual,
            double tolerance) {

        if (expected.length != actual.length) {
            throw new AssertionError(
                    "Vector lengths differ: expected "
                            + expected.length + " but got " + actual.length);
        }

        for (int index = 0; index < expected.length; index++) {
            if (!Double.isFinite(actual[index])
                    || Math.abs(expected[index] - actual[index]) > tolerance) {
                throw new AssertionError(
                        "Expected " + expected[index]
                                + " but got " + actual[index]
                                + " at index " + index);
            }
        }
    }

    private static void assertMatrixClose(double[][] expected, double[][] actual) {
        for (int row = 0; row < expected.length; row++) {
            for (int column = 0; column < expected[row].length; column++) {
                if (!Double.isFinite(actual[row][column])
                        || Math.abs(expected[row][column] - actual[row][column]) > TOLERANCE) {
                    throw new AssertionError("Expected " + expected[row][column]
                            + " but got " + actual[row][column]);
                }
            }
        }
    }

    private static void assertMatrixRelativeClose(
            double[][] expected, double[][] actual, double relativeTolerance) {

        double scale = 0.0;
        for (double[] row : expected) {
            for (double value : row) {
                scale = Math.max(scale, Math.abs(value));
            }
        }
        double tolerance = Math.max(Double.MIN_NORMAL, scale * relativeTolerance);

        for (int row = 0; row < expected.length; row++) {
            for (int column = 0; column < expected[row].length; column++) {
                if (!Double.isFinite(actual[row][column])
                        || Math.abs(expected[row][column] - actual[row][column]) > tolerance) {
                    throw new AssertionError("Expected " + expected[row][column]
                            + " but got " + actual[row][column]
                            + " with tolerance " + tolerance);
                }
            }
        }
    }
}

