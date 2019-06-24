package net.maizegenetics.analysis.modelfitter

import krangl.*
import krangl.experimental.oneHot
import net.maizegenetics.analysis.data.ExportPlugin
import net.maizegenetics.dna.snp.GenotypeTableUtils
import net.maizegenetics.dna.snp.ImportUtils
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition
import net.maizegenetics.stats.linearmodels.LinearModelUtils
import net.maizegenetics.matrixalgebra.Matrix.*
import net.maizegenetics.phenotype.PhenotypeBuilder
import net.maizegenetics.phenotype.PhenotypeUtils
import net.maizegenetics.tassel.TasselLogging
import net.maizegenetics.util.LoggingUtils
import kotlin.math.pow

fun main() {
    LoggingUtils.setupDebugLogging()
    /*
    val geno = ImportUtils.readFromHapmap("/Users/SAMUELFERNANDES/Documents/UIUC2/tassel-5-source/data_for_tests/genotype.hmp.txt", null, true)
    val pheno = PhenotypeBuilder().fromFile("/Users/SAMUELFERNANDES/Documents/UIUC2/tassel-5-source/data_for_tests/pheno.txt").build()[0]
    val Y = DoubleMatrixFactory.DEFAULT.make(pheno.numberOfObservations())
    val X = geno.genotypeMatrix()
    */

    val data = DataFrame.readTSV("/Users/SAMUELFERNANDES/Documents/UIUC2/tassel-5-source/data_for_tests/pheno.txt")
    val design = DataFrame.readCSV("/Users/SAMUELFERNANDES/Documents/UIUC2/tassel-5-source/data_for_tests/Design.hmp.csv")

    val Y = EJMLDoubleMatrix(data.remove("<Trait>").toDoubleMatrix()).transpose()

    //var test = design.names.groupingBy {
    //    it.substring(0, it.indexOf('_'))
    //}

    println(design.names)

    val xMu = EJMLDoubleMatrix(design.select("intercept").toDoubleMatrix()).transpose()
    val xMuSnp1 = EJMLDoubleMatrix(design.select("intercept", "snp1_GG", "snp1_AA").toDoubleMatrix()).transpose()
    val xMuSnp1Snp2 = EJMLDoubleMatrix(design.select("intercept", "snp1_GG", "snp1_AA", "snp2_GG", "snp2_AA").toDoubleMatrix()).transpose()
    val xMuSnp1Snp2Snp3 = EJMLDoubleMatrix(design.select("intercept", "snp1_GG", "snp1_AA", "snp2_GG", "snp2_AA", "snp3_AA", "snp3_CC").toDoubleMatrix()).transpose()
    val xMuSnp1Snp2Snp3Snp4 = EJMLDoubleMatrix(design.select("intercept", "snp1_GG", "snp1_AA", "snp2_GG", "snp2_AA", "snp3_AA", "snp3_CC", "snp4_AA", "snp4_GG").toDoubleMatrix()).transpose()
    val xMuSnp1Snp2Snp3Snp4Snp5 = EJMLDoubleMatrix(design.select("intercept", "snp1_GG", "snp1_AA", "snp2_GG", "snp2_AA", "snp3_AA", "snp3_CC", "snp4_AA", "snp4_GG", "snp5_AA", "snp5_GG", "snp5_AG").toDoubleMatrix()).transpose()
    val xFull = EJMLDoubleMatrix(design.toDoubleMatrix()).transpose()

    val pvalueI = manovaPvalue(Y, xMu)
    val pvalue1 = manovaPvalue(Y, xMuSnp1, xMu)
    val pvalue2 = manovaPvalue(Y, xMuSnp1Snp2, xMuSnp1)
    val pvalue3 = manovaPvalue(Y, xMuSnp1Snp2Snp3, xMuSnp1Snp2)
    val pvalue4 = manovaPvalue(Y, xMuSnp1Snp2Snp3Snp4, xMuSnp1Snp2Snp3)
    val pvalue5 = manovaPvalue(Y, xMuSnp1Snp2Snp3Snp4Snp5, xMuSnp1Snp2Snp3Snp4)
    val pvalue6 = manovaPvalue(Y, xFull, xMuSnp1Snp2Snp3Snp4Snp5)

    println("$pvalueI \n $pvalue1 \n $pvalue2 \n $pvalue3 \n $pvalue4 \n $pvalue5 \n $pvalue6")
}

data class BetaValue(val B: DoubleMatrix, val H: DoubleMatrix)
fun calcBeta(Y: DoubleMatrix, X: DoubleMatrix): BetaValue {
    val XtX: DoubleMatrix = X.crossproduct()
    val XtY: DoubleMatrix = X.crossproduct(Y)
    val XtXinv: DoubleMatrix = XtX.generalizedInverse()
    val B: DoubleMatrix = XtXinv.mult(XtY)
    val H: DoubleMatrix = B.transpose().mult(XtY)
    return BetaValue(B, H)
}

// Y = Phenotypic data
// xF = Full model
// xR = Reduced model
// returns Wilk's lambda pvalue
fun manovaPvalue(Y: DoubleMatrix, xF: DoubleMatrix, xR: DoubleMatrix): Double {
    //total SQ
    val YtY: DoubleMatrix = Y.crossproduct()
    //number of variables
    val p = Y.numberOfColumns().toDouble()
    //full model
    val hF = calcBeta(Y, xF).H
    //Residual (Error) matrix
    val E: DoubleMatrix = YtY.minus(hF)
    //reduced model
    val hR = calcBeta(Y, xR).H
    //adjusted hypothesis matrix
    val hA: DoubleMatrix = hF.minus(hR)
    //Wilks lambda
    val eigen: EigenvalueDecomposition = E.inverse().mult(hA).eigenvalueDecomposition
    var lambda = 1.0
    for (i in 0..(eigen.eigenvalues.size - 1)) {
        lambda *= 1 / (1 + eigen.eigenvalues[i])
    }
    //Degree of Freedom
    val df = xF.columnRank().toDouble() - xR.columnRank().toDouble()//data[data.names[0]].asStrings().distinctBy {it.hashCode()}.size.toDouble()
    val t: Double
    if ((p.pow(2) + df.pow(2) - 5) > 0) {
        t = Math.sqrt((p.pow(2) * df.pow(2) - 4) / (p.pow(2) + df.pow(2) - 5))
    } else {
        t = 1.0
    }

    val ve = Y.numberOfRows().toDouble() - xF.columnRank().toDouble()

    val r = ve - (p - df + 1) / 2
    val f = (p * df - 2) / 4
    val fCalc = ((1 - lambda.pow(1 / t)) / (lambda.pow(1 / t))) * ((r * t - 2 * f) / (p * df))
    val num: Double = p * df
    val den: Double = (r * t) - (2 * f)
    val pvalue: Double = LinearModelUtils.Ftest(fCalc, num, den)
    return pvalue
    //return listOf(ve, den, pvalue, lambda)
}

// Y = Phenotypic data
// xI = intercept
// returns Wilk's lambda pvalue
fun manovaPvalue(Y: DoubleMatrix, xI: DoubleMatrix): Double {
    //total SQ
    val YtY: DoubleMatrix = Y.crossproduct()
    //number of variables
    val p = Y.numberOfColumns().toDouble()
    //full model
    val hI = calcBeta(Y, xI).H
    //Residual (Error) matrix
    val E: DoubleMatrix = YtY.minus(hI)
    //Wilks lambda
    val eigen: EigenvalueDecomposition = E.inverse().mult(hI).eigenvalueDecomposition
    var lambda = 1.0
    for (i in 0..(eigen.eigenvalues.size - 1)) {
        lambda *= 1 / (1 + eigen.eigenvalues[i])
    }
    //Degree of Freedom
    val df = xI.columnRank().toDouble()//data[data.names[0]].asStrings().distinctBy {it.hashCode()}.size.toDouble()
    val t: Double
    if ((p.pow(2) + df.pow(2) - 5) > 0) {
        t = Math.sqrt((p.pow(2) * df.pow(2) - 4) / (p.pow(2) + df.pow(2) - 5))
    } else {
        t = 1.0
    }

    val ve = Y.numberOfRows().toDouble() - xI.columnRank().toDouble()

    val r = ve - (p - df + 1) / 2
    val f = (p * df - 2) / 4
    val fCalc = ((1 - lambda.pow(1 / t)) / (lambda.pow(1 / t))) * ((r * t - 2 * f) / (p * df))
    val num: Double = p * df
    val den: Double = (r * t) - (2 * f)
    val pvalue: Double = LinearModelUtils.Ftest(fCalc, num, den)

    return pvalue
    //return listOf(ve, den, pvalue, lambda)
}