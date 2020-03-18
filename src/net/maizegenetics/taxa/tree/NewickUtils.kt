@file:JvmName("NewickUtils")

package net.maizegenetics.taxa.tree

import org.apache.log4j.Logger
import java.io.BufferedWriter
import java.io.File
import java.lang.Double
import java.util.*
import kotlin.collections.HashMap

/**
 * These utilities are related to the Newick Tree Format.
 * http://evolution.genetics.washington.edu/phylip/newicktree.html
 */

private val myLogger = Logger.getLogger("net.maizegenetics.taxa.tree.NewickUtils")

/**
 * Creates a Tree from the given newick formatted file
 */
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

/**
 * Creates a node in the tree recursively creates the children for that node.
 */
private fun makeNode(newick: String): Node {

    val nameBranchLength = newick.substringAfterLast(')')
    val name = nameBranchLength.substringBefore(':').replace("'", "").let {
        try {
            Double.parseDouble(it)
            ""
        } catch (ne: NumberFormatException) {
            it.replace("_", "")
        }
    }
    val branchLength = if (nameBranchLength.contains(':')) nameBranchLength.substringAfter(':').toDouble() else 0.0

    val result = SimpleNode(name, branchLength)

    if (nameBranchLength != newick) {

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

/**
 * Writes give tree to Newick formatted file.
 * http://evolution.genetics.washington.edu/phylip/newicktree.html
 */
fun write(filename: String, tree: Tree, includeBranchLengths: Boolean = true) {

    try {
        File(filename).bufferedWriter().use { writer ->
            write(tree.root, writer, includeBranchLengths)
            writer.append(";\n")
        }
    } catch (e: Exception) {
        myLogger.debug(e.message, e)
        throw IllegalStateException("NewickUtils: write: problem writing: $filename.\n${e.message}")
    }

}

private fun write(node: Node, writer: BufferedWriter, includeBranchLengths: Boolean) {

    if (!node.isLeaf) {
        writer.append("(")
        for (i in 0 until node.childCount) {
            if (i != 0) writer.append(",")
            write(node.getChild(i), writer, includeBranchLengths)
        }
        writer.append(")")
    }

    node.identifier?.name?.let {
        if (it.isNotEmpty()) {
            writer.append("'")
            writer.append(it)
            writer.append("'")
        }
    }

    if (includeBranchLengths && node.branchLength != 0.0) {
        writer.append(":")
        writer.append(node.branchLength.toString())
    }

}

private const val MERGE_ROOT_NODE = "MERGE_ROOT_NODE"

fun mergeTrees(trees: List<Tree>): Tree {

    if (trees.size < 2) {
        throw IllegalArgumentException("NewickUtils: mergeTrees: must supply at least 2 trees.")
    }

    val nameToNode = HashMap<String, Node>()
    val nodeToNode = IdentityHashMap<Node, Node>()

    lateinit var rootNode: Node

    trees.forEach { tree ->

        val nodes = tree.nodes()

        nodes.forEach { node ->

            val nodeName = node.identifier.name
            val existingNode = nameToNode[nodeName]
            val existingParentNode = nodeToNode[node.parent]
            val existingParentName = existingParentNode?.identifier?.name
            val parentName = node.parent?.identifier?.name

            if (nodeToNode.isEmpty()) {
                rootNode = SimpleNode(nodeName, node.branchLength)
                nodeToNode[node] = rootNode
                if (nodeName.isNotEmpty()) nameToNode[nodeName] = rootNode!!
            } else if (existingNode != null) {
                if (existingParentName != null && parentName != null && existingParentName != parentName) {
                    throw IllegalArgumentException("NewickUtils: mergeTrees: node: $nodeName has different parents ($parentName, ${existingParentName}) in different trees.")
                }
                if (existingNode.branchLength != node.branchLength) {
                    myLogger.warn("mergeTrees: nodes named: $nodeName in different trees have different branch lengths.")
                }
                nodeToNode[node] = existingNode
            } else if (existingParentNode == null) {
                val mergeRoot = if (nameToNode[MERGE_ROOT_NODE] != null) {
                    nameToNode[MERGE_ROOT_NODE]!!
                } else {
                    val temp = SimpleNode("", 0.0)
                    nameToNode[MERGE_ROOT_NODE] = temp
                    temp
                }
                if (rootNode != mergeRoot) mergeRoot.addChild(rootNode)
                val newChild = SimpleNode(nodeName, node.branchLength)
                nodeToNode[node] = newChild
                if (nodeName.isNotEmpty()) nameToNode[nodeName] = newChild
                mergeRoot.addChild(newChild)

                rootNode = mergeRoot
            } else {
                val newChild = SimpleNode(nodeName, node.branchLength)
                nodeToNode[node] = newChild
                if (nodeName.isNotEmpty()) nameToNode[nodeName] = newChild
                existingParentNode.addChild(newChild)
            }

        }

    }

    return SimpleTree(rootNode)

}

fun Tree.nodes(): List<Node> {

    val result = mutableListOf<Node>()
    result.add(this.root)
    addChildren(result, this.root)
    return result

}

private fun addChildren(nodes: MutableList<Node>, node: Node) {

    for (i in 0 until node.childCount) {
        nodes.add(node.getChild(i))
        addChildren(nodes, node.getChild(i))
    }

}