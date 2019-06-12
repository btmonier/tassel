package net.maizegenetics.analysis.modelfitter

import krangl.*
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition
import net.maizegenetics.stats.linearmodels.LinearModelUtils
import net.maizegenetics.matrixalgebra.Matrix.*
import kotlin.math.pow


fun main() {
    var data = DataFrame.readCSV("/Users/pjbrown/Desktop/UIUC2/tassel-5-source/data_for_tests/fatorial_multi.csv")
    val design = DataFrame.readCSV("/Users/pjbrown/Desktop/UIUC2/tassel-5-source/data_for_tests/design_fatorial.csv")
    data = data.addColumn("interaction") { it["gender"] + "" + it["level"] }
    val Y = EJMLDoubleMatrix(data.select("k1", "k2", "k3").toDoubleMatrix()).transpose()
    val X = EJMLDoubleMatrix(design.toDoubleMatrix()).transpose()
    println(data.names)
    val xMu = EJMLDoubleMatrix(design.select("intercept").toDoubleMatrix()).transpose()
    val xMuAlpha = EJMLDoubleMatrix(design.select("intercept", "j1", "j2").toDoubleMatrix()).transpose()
    val xMuBeta = EJMLDoubleMatrix(design.select("intercept", "Y1", "Y2", "Y3", "Y4").toDoubleMatrix()).transpose()
    val xMuAlphaBeta = EJMLDoubleMatrix(design.remove("j1Y1", "j2Y1", "j1Y2", "j2Y2", "j1Y3", "j2Y3", "j1Y4", "j2Y4").toDoubleMatrix()).transpose()
    val xTestAlpha = EJMLDoubleMatrix(design.remove("j1", "j2").toDoubleMatrix()).transpose()

    val pvalue1: Double = manovaPvalue(Y, xMuAlpha, xMu)
    val pvalue2: Double = manovaPvalue(Y, xMuBeta, xMu)
    val pvalue3: Double = manovaPvalue(Y, xMuAlphaBeta, xMuBeta)
    val pvalue4: Double = manovaPvalue(Y, xMuAlphaBeta, xMuAlpha)
    val pvalue5: Double = manovaPvalue(Y, X, xMuAlphaBeta)
    //val pvalueT: Double = manovaPvalue(Y, X, xTestAlpha)

    println("$pvalue1 \n $pvalue2 \n $pvalue3 \n $pvalue4 \n $pvalue5 ")
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

//Y = Phenotypic data
//xR = Reduced model
//xF = Full model
// returns Wilk's lambda pvalue
fun manovaPvalue(Y: DoubleMatrix, xF: DoubleMatrix, xR: DoubleMatrix): Double{
    //total SQ
    val YtY: DoubleMatrix = Y.crossproduct()
    //number of variables
    var p = Y.numberOfColumns().toDouble()
    //full model
    val (bF, hF) = calcBeta(Y, xF)
    //reduced model
    val (bR, hR) = calcBeta(Y, xR)
    //adjusted hypotesis matrix
    val hA : DoubleMatrix = hF.minus(hR)
    //Residual (Error) matrix
    val E : DoubleMatrix = YtY.minus(hF)
    //Wilks lambda
    val eigen : EigenvalueDecomposition =  E.inverse().mult(hA).eigenvalueDecomposition
    var lambda  = 1.0
    for(i in 0..(eigen.eigenvalues.size-1) ){
        lambda *= 1 / (1 + eigen.eigenvalues[i])
    }
    //Degree of Freedom
    val df = xF.columnRank().toDouble() - xR.columnRank().toDouble()//data[data.names[0]].asStrings().distinctBy {it.hashCode()}.size.toDouble()
    val t: Double
    if((p.pow(2) + df.pow(2) - 5) > 0){
        t = Math.sqrt((p.pow(2)*df.pow(2) - 4)/(p.pow(2) + df.pow(2) -5))
    }else{t = 1.0}

    val ve =Y.numberOfRows().toDouble() - df - 1

    val r = ve - (p-df+1)/2
    val f= (p*df-2)/4
    val fCalc = ((1-lambda.pow(1/t))/(lambda.pow(1/t)))*((r*t - 2*f)/(p*df))
    val pvalue : Double = LinearModelUtils.Ftest(fCalc, p*df, (r*t)-(2*f))

    return pvalue
}