public class QRDecompositionTest {
    private static final double TOLERANCE = 1e-9;

    public static void main(String[] args) {
        testQrRebuildsOriginalMatrix();
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
}
