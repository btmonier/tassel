package net.maizegenetics.matrixalgebra.Matrix

import krangl.*
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition
import net.maizegenetics.stats.linearmodels.LinearModelUtils
import kotlin.math.pow
//import org.apache.commons.math3.stat.inference.OneWayAnova



fun main() {
    val data = DataFrame.readCSV("/Users/SAMUELFERNANDES/Documents/UIUC2/tassel-5-source/data_for_tests/multi.csv")
    val design = DataFrame.readCSV("/Users/SAMUELFERNANDES/Documents/UIUC2/tassel-5-source/data_for_tests/design.csv")
    //var X : Array<DoubleArray> = design.toDoubleMatrix()
    val X = EJMLDoubleMatrix(design.toDoubleMatrix()).transpose()
    //intercept only
    val intercept = EJMLDoubleMatrix(design.select("intercept").toDoubleMatrix()).transpose()
    val pheno = EJMLDoubleMatrix(data.remove("taxa").toDoubleMatrix()).transpose()

    val YtY: DoubleMatrix = pheno.crossproduct()
    val (B, R, E) = calcBeta(pheno, X)
    val (B2, R2, E2) = calcBeta(pheno, intercept)

    val H : DoubleMatrix = R.minus(R2)
    val Total : DoubleMatrix = YtY.minus(R2)

    val eigen : EigenvalueDecomposition =  E.inverse().mult(H).eigenvalueDecomposition

    var lambda  = 1.0
    for(i in 0..(eigen.eigenvalues.size-1) ){
        lambda *= 1 / (1 + eigen.eigenvalues[i])
    }

    var p = pheno.numberOfColumns()
    val vh = 3
    val ve = pheno.numberOfRows() - (vh + 1)
    val t : Double

    if((p.toDouble().pow(2) + vh.toDouble().pow(2) - 5 ) > 0) {
        t = Math.sqrt((p.toDouble().pow(2) * vh.toDouble().pow(2) - 4) / (p.toDouble().pow(2) + vh.toDouble().pow(2) - 5))
    } else {
        t = 1.0
    }

    val r : Int = ve - (p - vh + 1)
    val f : Int = (p * vh - 2)/4

    val Fc : Double = ((1 - (lambda.pow(1/t)))/(lambda.pow(1/t))) * ((r*t - 2*f)/(p*vh))

    val pvalue : Double = LinearModelUtils.Ftest(Fc, p*vh.toDouble(), (r*t)-(2*f))

    println(pvalue)

 /*
    val n: Int = eigen.eigenvalues.size
    //fun f(x: Double) = (1..100).map { i -> 3 * x.pow(2) + (2*i) }.sum()
    val list = arrayListOf<Int>(1,2,3,4)
    println((1..n).fold(list, Long::times))

 //Ftest
 inline fun <T, R> Array<out T>.fold( initial: R, operation: (acc: R, T) -> R ): R
*/
     }

data class BetaValue(val B: DoubleMatrix, val R: DoubleMatrix, val E: DoubleMatrix)

fun calcBeta(Y: DoubleMatrix, X: DoubleMatrix): BetaValue {
    val YtY: DoubleMatrix = Y.crossproduct()
    val XtX: DoubleMatrix = X.crossproduct()
    val XtY: DoubleMatrix = X.crossproduct(Y)
    val XtXinv: DoubleMatrix = XtX.generalizedInverse()
    val B: DoubleMatrix = XtXinv.mult(XtY)
    val R: DoubleMatrix = B.transpose().mult(XtY)
    val E: DoubleMatrix = YtY.minus(R)
    return BetaValue(B, R, E)
}



/**
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


import java.util.ArrayList

import org.apache.commons.math3.exception.MathIllegalArgumentException
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import org.apache.commons.math3.stat.inference.OneWayAnova
import org.junit.Assert
import org.junit.Test



 * Test cases for the OneWayAnovaImpl class.
 *


class OneWayAnovaTest {

    protected var testStatistic = OneWayAnova()

    private val emptyArray = doubleArrayOf()

    private val classA = doubleArrayOf(93.0, 103.0, 95.0, 101.0, 91.0, 105.0, 96.0, 94.0, 101.0)
    private val classB = doubleArrayOf(99.0, 92.0, 102.0, 100.0, 102.0, 89.0)
    private val classC = doubleArrayOf(110.0, 115.0, 111.0, 117.0, 128.0, 117.0)

    @Test
    fun testAnovaFValue() {
        // Target comparison values computed using R version 2.6.0 (Linux version)
        val threeClasses: List<ArrayList> = ArrayList()
        threeClasses.add(classA)
        threeClasses.add(classB)
        threeClasses.add(classC)

        Assert.assertEquals("ANOVA F-value", 24.67361709460624,
                testStatistic.anovaFValue(threeClasses), 1E-12)

        val twoClasses = ArrayList()
        twoClasses.add(classA)
        twoClasses.add(classB)

        Assert.assertEquals("ANOVA F-value", 0.0150579150579,
                testStatistic.anovaFValue(twoClasses), 1E-12)

        val emptyContents = ArrayList()
        emptyContents.add(emptyArray)
        emptyContents.add(classC)
        try {
            testStatistic.anovaFValue(emptyContents)
            Assert.fail("empty array for key classX, MathIllegalArgumentException expected")
        } catch (ex: MathIllegalArgumentException) {
            // expected
        }

        val tooFew = ArrayList()
        tooFew.add(classA)
        try {
            testStatistic.anovaFValue(tooFew)
            Assert.fail("less than two classes, MathIllegalArgumentException expected")
        } catch (ex: MathIllegalArgumentException) {
            // expected
        }

    }


    @Test
    fun testAnovaPValue() {
        // Target comparison values computed using R version 2.6.0 (Linux version)
        val threeClasses = ArrayList()
        threeClasses.add(classA)
        threeClasses.add(classB)
        threeClasses.add(classC)

        Assert.assertEquals("ANOVA P-value", 6.959446E-06,
                testStatistic.anovaPValue(threeClasses), 1E-12)

        val twoClasses = ArrayList()
        twoClasses.add(classA)
        twoClasses.add(classB)

        Assert.assertEquals("ANOVA P-value", 0.904212960464,
                testStatistic.anovaPValue(twoClasses), 1E-12)

    }

    @Test
    fun testAnovaPValueSummaryStatistics() {
        // Target comparison values computed using R version 2.6.0 (Linux version)
        val threeClasses = ArrayList()
        val statsA = SummaryStatistics()
        for (a in classA) {
            statsA.addValue(a)
        }
        threeClasses.add(statsA)
        val statsB = SummaryStatistics()
        for (b in classB) {
            statsB.addValue(b)
        }
        threeClasses.add(statsB)
        val statsC = SummaryStatistics()
        for (c in classC) {
            statsC.addValue(c)
        }
        threeClasses.add(statsC)

        Assert.assertEquals("ANOVA P-value", 6.959446E-06,
                testStatistic.anovaPValue(threeClasses, true), 1E-12)

        val twoClasses = ArrayList()
        twoClasses.add(statsA)
        twoClasses.add(statsB)

        Assert.assertEquals("ANOVA P-value", 0.904212960464,
                testStatistic.anovaPValue(twoClasses, false), 1E-12)

    }

    @Test
    fun testAnovaTest() {
        // Target comparison values computed using R version 2.3.1 (Linux version)
        val threeClasses = ArrayList()
        threeClasses.add(classA)
        threeClasses.add(classB)
        threeClasses.add(classC)

        Assert.assertTrue("ANOVA Test P<0.01", testStatistic.anovaTest(threeClasses, 0.01))

        val twoClasses = ArrayList()
        twoClasses.add(classA)
        twoClasses.add(classB)

        Assert.assertFalse("ANOVA Test P>0.01", testStatistic.anovaTest(twoClasses, 0.01))
    }

}
 */