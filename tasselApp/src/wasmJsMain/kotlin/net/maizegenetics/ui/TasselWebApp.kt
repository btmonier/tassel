package net.maizegenetics.ui

import androidx.compose.animation.AnimatedVisibility
import androidx.compose.foundation.border
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.unit.dp

@Composable
fun TasselWebApp() {
    MaterialTheme {

        Row(
            modifier = Modifier
                .fillMaxSize()
                .padding(16.dp),
            horizontalArrangement = Arrangement.Center,
            verticalAlignment = Alignment.Top
        ) {
            leftColumn()
            Spacer(modifier = Modifier.width(10.dp))
            exampleColumn()
        }

    }
}

@Composable
private fun leftColumn() {
    Column(
        modifier = Modifier
            .safeContentPadding()
            .fillMaxSize(0.25f)
            .border(2.dp, Color.Red, shape = RoundedCornerShape(4.dp))
            .padding(8.dp),
        horizontalAlignment = Alignment.CenterHorizontally,
    ) {

        Row(
            modifier = Modifier
                .fillMaxWidth()
                .border(2.dp, Color.Red, shape = RoundedCornerShape(4.dp))
                .padding(8.dp)
                .weight(1f),
            horizontalArrangement = Arrangement.SpaceBetween
        ) {
            Text("This is a Tassel Web App")
        }

        var text by remember { mutableStateOf("Initial value") }

        Row(
            modifier = Modifier
                .fillMaxWidth()
                .border(2.dp, Color.Red, shape = RoundedCornerShape(4.dp))
                .padding(8.dp)
                .weight(1f),
            horizontalArrangement = Arrangement.SpaceBetween
        ) {
            TextField(
                value = text,
                onValueChange = { /* no-op: ignore user edits */ },
                readOnly = true,
                label = { Text("Dataset Description") },
                modifier = Modifier
                    .fillMaxWidth()
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