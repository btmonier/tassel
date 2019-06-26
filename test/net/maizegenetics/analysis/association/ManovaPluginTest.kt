package net.maizegenetics.analysis.association

import org.junit.Test

class ManovaPluginTest {

    @Test
    fun imputeNTest() {
        val genotype = arrayOf("A","A","A","A","A","A","T","T","N","N","N","N","N","N","N","N")
        println("${genotype.joinToString(",")}")

        println("${ManovaPlugin(null, false).imputeNsInGenotype(genotype).joinToString(",")}")
    }
}