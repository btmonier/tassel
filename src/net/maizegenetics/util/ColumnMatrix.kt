package net.maizegenetics.util

/**
 * @author Terry Casstevens
 * Created December 06, 2018
 */

class ColumnMatrix private constructor(private val matrix: SuperByteMatrix) {

    val numColumns = matrix.numRows
    val numRows = matrix.numColumns

    fun column(index: Int) = matrix.getAllColumns(index)

    class Builder(numRows: Int, numColumns: Int) {

        private val matrix = SuperByteMatrixBuilder.getInstance(numColumns, numRows)

        fun set(row: Int, column: Int, value: Byte) = matrix.set(column, row, value)

        fun build() = ColumnMatrix(matrix)
    }

}

fun main(args: Array<String>) {
    var numRows = 5
    var numColumns = 10
    var builder = ColumnMatrix.Builder(numRows, numColumns)
    for (r in 0 until numRows) {
        for (c in 0 until numColumns) {
            builder.set(r, c, (r + c).toByte())
        }
    }
    val matrix = builder.build()
    for (c in 0 until numColumns) {
        val column = matrix.column(c)
        for (r in 0 until numRows) {
            print("${column[r]} ")
        }
        println()
    }
}