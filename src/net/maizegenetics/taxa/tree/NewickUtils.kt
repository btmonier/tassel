@file:JvmName("NewickUtils")

package net.maizegenetics.taxa.tree

import java.io.File

fun read(filename: String): Tree {

    val newick = File(filename).readLines().joinToString(separator = "") { it }

    if (newick.count { it == '(' } != newick.count { it == ')' }) {
        throw IllegalArgumentException("NewickUtils: read: $filename: number of open parenthesis doesn't match number of close parenthesis.")
    }

    if (newick.last() != ';') {
        throw IllegalArgumentException("NewickUtils: read: $filename: doesn't end with semicolon.")
    }

    return SimpleTree(makeNode(newick.substring(0, newick.length - 1)))

}

private fun makeNode(newick: String): Node {

    val nameWeight = newick.substringAfterLast(')')
    val name = nameWeight.substringBefore(':').replace("'", "")
    val weight = if (nameWeight.contains(':')) nameWeight.substringAfter(':').toDouble() else 0.0

    val result = SimpleNode(name, weight)

    if (nameWeight != newick) {

        var parenthesisCount = 0
        var currentStr = StringBuilder()
        newick.substringAfter('(').substringBeforeLast(')').forEach {
            when (it) {
                '(' -> {
                    parenthesisCount++
                    currentStr.append(it)
                }
                ')' -> {
                    parenthesisCount--
                    currentStr.append(it)
                }
                ',' -> {
                    if (parenthesisCount == 0) {
                        result.addChild(makeNode(currentStr.toString()))
                        currentStr = StringBuilder()
                    } else {
                        currentStr.append(it)
                    }
                }
                else -> currentStr.append(it)
            }
        }
        if (currentStr.isNotEmpty()) result.addChild(makeNode(currentStr.toString()))

    }

    return result

}


fun main() {
    val filename = "/Users/tmc46/projects/gerp_pipeline/phytozome_12.nwk"
    read(filename)
}