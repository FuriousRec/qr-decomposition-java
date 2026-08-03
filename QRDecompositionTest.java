public class QRDecompositionTest {
    private static final double TOLERANCE = 1e-9;

    public static void main(String[] args) {
        testQrRebuildsOriginalMatrix();
        System.out.println("Test passed.");
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
                if (Math.abs(expected[row][column] - actual[row][column]) > TOLERANCE) {
                    throw new AssertionError("Expected " + expected[row][column]
                            + " but got " + actual[row][column]);
                }
            }
        }
    }
}
