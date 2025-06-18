package net.maizegenetics.ui

import androidx.compose.animation.AnimatedVisibility
import androidx.compose.foundation.Image
import androidx.compose.foundation.border
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.filled.ChevronRight
import androidx.compose.material.icons.filled.ExpandMore
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.unit.Dp
import androidx.compose.ui.unit.dp
import org.jetbrains.compose.resources.painterResource
import tassel.tasselapp.generated.resources.Res
import tassel.tasselapp.generated.resources.genotype_viewer_example

//import net.maizegenetics.dna.snp.ImportUtils

@Composable
fun TasselWebApp() {

    //val genotype = ImportUtils.readFromVCF("/Users/tmc46/git/tassel-5-standalone/TASSELTutorialData/data/mdp_genotype.hmp.txt", null, false, true)

    MaterialTheme {

        Row(
            modifier = Modifier
                .fillMaxSize()
                .fillMaxHeight()
                .fillMaxWidth()
                .padding(10.dp)
                .border(2.dp, Color.Gray, shape = RoundedCornerShape(10.dp)),
            horizontalArrangement = Arrangement.Center,
            verticalAlignment = Alignment.Top
        ) {
            leftColumn()
            Spacer(modifier = Modifier.width(10.dp))
            // exampleColumn()
            ResourceImageExample()
            // rightColumn()
        }

    }

}

@Composable
private fun leftColumn() {
    Column(
        modifier = Modifier
            .safeContentPadding()
            .fillMaxWidth(0.25f)
            .fillMaxHeight()
            //.border(2.dp, Color.Red, shape = RoundedCornerShape(7.dp))
            .padding(8.dp),
        horizontalAlignment = Alignment.CenterHorizontally,
    ) {
        datasetTree(modifier = Modifier.weight(0.45f))
        datasetInfo(modifier = Modifier.weight(0.4f))
        progressBar(modifier = Modifier.weight(0.15f))
    }
}

@Composable
private fun datasetTree(modifier: Modifier = Modifier) {
    Row(
        modifier = Modifier
            .fillMaxWidth()
            .fillMaxHeight(0.45f)
            //.border(2.dp, Color.Red, shape = RoundedCornerShape(4.dp))
            .padding(8.dp)
    ) {
        Column(
            modifier = Modifier
                .fillMaxWidth()
                .fillMaxHeight(),
            horizontalAlignment = Alignment.CenterHorizontally,
        ) {
            Text("Datasets")
            Spacer(modifier = Modifier.height(5.dp))
            DatasetTreeExample()
        }
    }
}

@Composable
private fun datasetInfo(modifier: Modifier = Modifier) {

    var text by remember {
        mutableStateOf(
            "Number of taxa: 281\n" +
                    "Number of sites: 3093\n" +
                    "\n" +
                    "Chromosomes...\n" +
                    "\n" +
                    "1: 540 sites:\n" +
                    "0 (157104) - 539 (299170077)\n" +
                    "\n" +
                    "2: 393 sites:\n" +
                    "540 (736367) - 932 (234574991)\n" +
                    "\n" +
                    "3: 355 sites:\n" +
                    "933 (1240310) - 1287 (229544509)\n" +
                    "\n" +
                    "4: 319 sites:\n" +
                    "1288 (139753) - 1606 (245131801)\n" +
                    "\n" +
                    "5: 357 sites:\n" +
                    "1607 (656148) - 1963 (216431558)\n" +
                    "\n" +
                    "6: 213 sites:\n" +
                    "1964 (2379148) - 2176 (167883450)\n" +
                    "\n" +
                    "7: 246 sites:\n" +
                    "2177 (729478) - 2422 (170346253)\n" +
                    "\n" +
                    "8: 256 sites:\n" +
                    "2423 (169137) - 2678 (172323795)\n" +
                    "\n" +
                    "9: 213 sites:\n" +
                    "2679 (3873116) - 2891 (151289948)\n" +
                    "\n" +
                    "10: 201 sites:\n" +
                    "2892 (838970) - 3092 (148907116)\n" +
                    "\n" +
                    "\n" +
                    "Nucleotide Codes\n" +
                    "(Derived from IUPAC)...\n" +
                    "A     A:A\n" +
                    "C     C:C\n" +
                    "G     G:G\n" +
                    "T     T:T\n" +
                    "R     A:G\n" +
                    "Y     C:T\n" +
                    "S     C:G\n" +
                    "W     A:T\n" +
                    "K     G:T\n" +
                    "M     A:C\n" +
                    "+     +:+ (insertion)\n" +
                    "0     +:-\n" +
                    "-     -:- (deletion)\n" +
                    "N     Unknown\n" +
                    "\n"
        )
    }

    Row(
        modifier = modifier
            .fillMaxWidth()
            .fillMaxHeight(0.4f)
            //.border(2.dp, Color.Red, shape = RoundedCornerShape(4.dp))
            .padding(8.dp)
    ) {
        Column(
            modifier = Modifier
                .fillMaxWidth()
                .fillMaxHeight(),
            horizontalAlignment = Alignment.CenterHorizontally,
        ) {
            Text("Dataset Description")
            Spacer(modifier = Modifier.height(10.dp))
            TextField(
                value = text,
                onValueChange = { /* no-op: ignore user edits */ },
                readOnly = true,
                //label = { Text("Dataset Description") },
                modifier = Modifier
                    .fillMaxWidth()
                    .fillMaxHeight()
            )
        }
    }

}

@Composable
private fun progressBar(modifier: Modifier = Modifier) {
    Row(
        modifier = Modifier
            .fillMaxWidth()
            .fillMaxHeight(0.15f)
            //.border(2.dp, Color.Red, shape = RoundedCornerShape(7.dp))
            .padding(8.dp)
    ) {
        Column(
            modifier = Modifier
                .fillMaxWidth()
                .fillMaxHeight(),
            horizontalAlignment = Alignment.CenterHorizontally,
        ) {
            Text("Progress")
            Spacer(modifier = Modifier.height(10.dp))
            LinearProgressIndicator(
                progress = {
                    0.25f
                },
                modifier = Modifier
                    .fillMaxWidth()
                    .height(10.dp)
            )
        }
    }
}

@Composable
private fun rightColumn() {
    TODO()
}

@Composable
fun exampleColumn() {
    var showContent by remember { mutableStateOf(false) }
    Column(
        modifier = Modifier
            .safeContentPadding()
            .fillMaxWidth(0.75f)
            .border(2.dp, Color.Red, shape = RoundedCornerShape(4.dp))
            .padding(8.dp),
        horizontalAlignment = Alignment.CenterHorizontally,
    ) {
        Button(onClick = { showContent = !showContent }) {
            Text("Click me!")
        }
        AnimatedVisibility(showContent) {
            val greeting = remember { "Hello Tassel" }
            Column(Modifier.fillMaxWidth(), horizontalAlignment = Alignment.CenterHorizontally) {
                //Image(painterResource(Res.drawable.compose_multiplatform), null)
                Text("Compose: $greeting")
                Text("Compose: $greeting")
            }
        }
    }
}

data class TreeNode<T>(
    val value: T,
    val children: List<TreeNode<T>> = emptyList()
)

@Composable
fun <T> TreeView(
    nodes: List<TreeNode<T>>,
    indent: Dp = 16.dp,
    nodeContent: @Composable (T) -> Unit
) {
    Column(
        modifier = Modifier
            .padding(8.dp)
            .fillMaxHeight()
        //.border(2.dp, Color.Red, shape = RoundedCornerShape(4.dp))
    ) {
        nodes.forEach { node ->
            TreeNodeView(node, 0, indent, nodeContent)
        }
    }
}

@Composable
private fun <T> TreeNodeView(
    node: TreeNode<T>,
    depth: Int,
    indent: Dp,
    nodeContent: @Composable (T) -> Unit
) {
    // track expand/collapse state per node
    var expanded by remember { mutableStateOf(false) }

    Column {
        Row(
            verticalAlignment = Alignment.CenterVertically,
            modifier = Modifier
                .fillMaxWidth()
                .padding(start = indent.times(depth), top = 4.dp, bottom = 4.dp)
                .clickable {
                    if (node.children.isNotEmpty()) expanded = !expanded
                }
        ) {
            // arrow icon for expand/collapse
            if (node.children.isNotEmpty()) {
                Icon(
                    imageVector = if (expanded) Icons.Default.ExpandMore else Icons.Default.ChevronRight,
                    contentDescription = if (expanded) "Collapse" else "Expand",
                    modifier = Modifier.size(16.dp)
                )
            } else {
                Spacer(modifier = Modifier.width(16.dp))
            }

            Spacer(Modifier.width(4.dp))

            // render the node's actual content
            nodeContent(node.value)
        }

        // if expanded, recurse into children
        if (expanded) {
            node.children.forEach { child ->
                TreeNodeView(child, depth + 1, indent, nodeContent)
            }
        }
    }
}

@Composable
private fun DatasetTreeExample() {
    val treeData = listOf(
        TreeNode(
            "Matrix", listOf(
                TreeNode("mdp_kinship"),
                // TreeNode("Dogs"),
                // TreeNode("Birds")
            )
        ),
        TreeNode("Numerical", listOf(
            TreeNode("mdp_phenotype")
        )),
        TreeNode("Sequence"),
        TreeNode("Results")
    )

    TreeView(nodes = treeData) { label ->
        Text(text = label, style = MaterialTheme.typography.bodyMedium)
    }
}

@Composable
fun ResourceImageExample() {

    Column(
        modifier = Modifier
            .safeContentPadding()
            .fillMaxWidth(0.75f)
            //.border(2.dp, Color.Red, shape = RoundedCornerShape(4.dp))
            .padding(8.dp),
        horizontalAlignment = Alignment.CenterHorizontally,
    ) {
        Image(painterResource(Res.drawable.genotype_viewer_example), null)
    }

}

