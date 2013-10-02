/*
 *  SuperByteMatrixSingle
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixSingle implements SuperByteMatrix {

    private final byte[] myData;
    private int myNumRows;
    private int myNumColumns;
    private final long myPrecompute1;

    SuperByteMatrixSingle(int rows, int columns) {

        myNumRows = rows;
        myNumColumns = columns;
        myPrecompute1 = (long) myNumColumns * (long) myNumRows - 1l;

        long numElements = (long) myNumRows * (long) myNumColumns;
        if (numElements > (long) (Integer.MAX_VALUE - 10)) {
            throw new IllegalArgumentException("SuperByteMatrixSingle: init: this number of rows: " + rows + "  and columns: " + columns + " is too large for SuperByteMatrixSingle.");
        }
        myData = new byte[(int) numElements];

    }

    @Override
    public void set(int row, int column, byte value) {
        myData[getIndex(row, column)] = value;
    }

    @Override
    public byte get(int row, int column) {
        return myData[getIndex(row, column)];
    }

    @Override
    public byte[] getAllColumns(int row) {

        if ((row < 0) || (row >= myNumRows)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingle: getAllColumns: row: " + row);
        }

        int start = getIndex(row, 0);
        byte[] result = new byte[myNumColumns];
        System.arraycopy(myData, start, result, 0, myNumColumns);
        return result;

    }

    @Override
    public byte[] getColumnRange(int row, int start, int end) {

        if ((row < 0) || (row >= myNumRows)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingle: getColumnRange: row: " + row);
        }

        if ((start < 0) || (start >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingle: getColumnRange: start: " + start);
        }

        if ((end < 0) || (end >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingle: getColumnRange: end: " + end);
        }

        if (end < start) {
            throw new IllegalArgumentException("SuperByteMatrixSingle: getColumnRange: end: " + end + " less than start: " + start);
        }

        int startIndex = getIndex(row, start);
        int numElements = end - start;
        byte[] result = new byte[numElements];
        System.arraycopy(myData, startIndex, result, 0, numElements);
        return result;

    }

    @Override
    public byte[] getAllRows(int column) {

        if ((column < 0) || (column >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingle: getAllRows: column: " + column);
        }

        byte[] result = new byte[myNumRows];
        int current = column;
        for (int i = 0; i < myNumRows; i++) {
            result[i] = myData[current];
            current += myNumColumns;
        }
        return result;

    }

    private int getIndex(int row, int column) {
        return row * myNumColumns + column;
    }

    @Override
    public int getNumRows() {
        return myNumRows;
    }

    @Override
    public int getNumColumns() {
        return myNumColumns;
    }

    @Override
    public boolean isColumnInnerLoop() {
        return true;
    }

    private int translateIndexForTranspose(long index) {
        return (int) ((index + (index % (long) myNumRows) * myPrecompute1) / (long) myNumRows);
    }

    public void transpose() {
        int numElements = myNumColumns * myNumRows;
        BitSet notVisited = new OpenBitSet(numElements);
        notVisited.set(0, numElements);
        int currentIndex = 1;
        byte temp;
        while (currentIndex != -1) {

            currentIndex = notVisited.nextSetBit(currentIndex);

            if (currentIndex != -1) {

                temp = myData[currentIndex];

                int srcIndex = translateIndexForTranspose(currentIndex);
                int destIndex = currentIndex;
                while (srcIndex != currentIndex) {
                    myData[destIndex] = myData[srcIndex];
                    notVisited.fastClear(destIndex);
                    destIndex = srcIndex;
                    srcIndex = translateIndexForTranspose(destIndex);
                }

                myData[destIndex] = temp;
                notVisited.fastClear(destIndex);

            }

        }

        int tempSize = myNumColumns;
        myNumColumns = myNumRows;
        myNumRows = tempSize;
    }

    @Override
    public void reorderRows(int[] newIndices) {

        if (newIndices.length != myNumRows) {
            throw new IllegalArgumentException("SuperByteMatrixSingle: reorderRows: index array size: " + newIndices.length + " doesn't equal num rows in matrix: " + myNumRows);
        }

        int[] tempIndices = new int[newIndices.length];
        System.arraycopy(newIndices, 0, tempIndices, 0, myNumRows);

        int currentRow = 0;
        byte[] temp = new byte[myNumColumns];

        while (currentRow < myNumRows) {

            while (currentRow < myNumRows) {
                if ((tempIndices[currentRow] == currentRow) || (tempIndices[currentRow] == -1)) {
                    tempIndices[currentRow] = -1;
                } else {
                    break;
                }
                currentRow++;
            }

            if (currentRow < myNumRows) {

                System.arraycopy(myData, getIndex(currentRow, 0), temp, 0, myNumColumns);

                int srcRow = tempIndices[currentRow];
                int destRow = currentRow;
                while (srcRow != currentRow) {
                    System.arraycopy(myData, getIndex(srcRow, 0), myData, getIndex(destRow, 0), myNumColumns);
                    tempIndices[destRow] = -1;
                    destRow = srcRow;
                    srcRow = tempIndices[destRow];
                }

                System.arraycopy(temp, 0, myData, getIndex(destRow, 0), myNumColumns);
                tempIndices[destRow] = -1;

            }

        }

    }

    @Override
    public void reorderColumns(int[] newIndices) {

        if (newIndices.length != myNumColumns) {
            throw new IllegalArgumentException("SuperByteMatrixSingle: reorderColumns: index array size: " + newIndices.length + " doesn't equal num columns in matrix: " + myNumColumns);
        }

        int[] tempIndices = new int[newIndices.length];
        System.arraycopy(newIndices, 0, tempIndices, 0, myNumColumns);

        int currentRow = 0;
        byte[] temp = new byte[myNumRows];

        while (currentRow < myNumColumns) {

            while (currentRow < myNumColumns) {
                if ((tempIndices[currentRow] == currentRow) || (tempIndices[currentRow] == -1)) {
                    tempIndices[currentRow] = -1;
                } else {
                    break;
                }
                currentRow++;
            }

            if (currentRow < myNumColumns) {

                for (int r = 0; r < myNumRows; r++) {
                    temp[r] = get(r, currentRow);
                }

                int srcColumn = tempIndices[currentRow];
                int destColumn = currentRow;
                while (srcColumn != currentRow) {
                    for (int r = 0; r < myNumRows; r++) {
                        set(r, destColumn, get(r, srcColumn));
                    }
                    tempIndices[destColumn] = -1;
                    destColumn = srcColumn;
                    srcColumn = tempIndices[destColumn];
                }

                for (int r = 0; r < myNumRows; r++) {
                    set(r, destColumn, temp[r]);
                }
                tempIndices[destColumn] = -1;

            }

        }

    }
}
