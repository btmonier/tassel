package net.maizegenetics.matrixalgebra.Matrix

import krangl.*
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition
import net.maizegenetics.stats.linearmodels.LinearModelUtils
import kotlin.math.pow
import java.util.HashSet
import java.util.ArrayList


//import org.apache.commons.math3.stat.inference.OneWayAnova

fun main() {
    var data = DataFrame.readCSV("/Users/pjbrown/Desktop/UIUC2/tassel-5-source/data_for_tests/fatorial_multi.csv")
    val design = DataFrame.readCSV("/Users/pjbrown/Desktop/UIUC2/tassel-5-source/data_for_tests/design_fatorial.csv")

    val X = EJMLDoubleMatrix(design.toDoubleMatrix()).transpose()

    data = data.addColumn("interaction") { it["gender"] + "" + it["level"] }

    val Y = EJMLDoubleMatrix(data.select("k1", "k2", "k3").toDoubleMatrix()).transpose()

    val X2 = EJMLDoubleMatrix(design.remove("j1Y1", "j2Y1", "j1Y2", "j2Y2", "j1Y3", "j2Y3", "j1Y4", "j2Y4").toDoubleMatrix()).transpose()
    val X3 = EJMLDoubleMatrix(design.select("intercept", "j1", "j2").toDoubleMatrix()).transpose()
    val X4 = EJMLDoubleMatrix(design.select("intercept", "Y1", "Y2", "Y3", "Y4").toDoubleMatrix()).transpose()
    val X5 = EJMLDoubleMatrix(design.select("intercept").toDoubleMatrix()).transpose()

    //mu, alpha, betha, delta
    val (B1, H1) = calcBeta(Y, X)
    //mu, alpha, betha
    val (B2, H2) = calcBeta(Y, X2)
    //mu, alpha
    val (B3, H3) = calcBeta(Y, X3)
    //mu, betha
    val (B4, H4) = calcBeta(Y, X4)
    //mu
    val (bInt, hInt) = calcBeta(Y, X5)

    val YtY: DoubleMatrix = Y.crossproduct()

    val hA : DoubleMatrix = H3.minus(hInt)
    val hB : DoubleMatrix = H4.minus(hInt)
    val hD : DoubleMatrix = H1.minus(H2)
    val Total : DoubleMatrix = YtY.minus(hA).minus(hB).minus(hD)

    val E : DoubleMatrix = YtY.minus(H1)

    val eigenInt : EigenvalueDecomposition =  E.inverse().mult(hInt).eigenvalueDecomposition
    var lambdaInt  = 1.0
    for(i in 0..(eigenInt.eigenvalues.size-1) ){
        lambdaInt *= 1 / (1 + eigenInt.eigenvalues[i])
    }

    val eigenA : EigenvalueDecomposition =  E.inverse().mult(hA).eigenvalueDecomposition
    var lambdaA  = 1.0
    for(i in 0..(eigenA.eigenvalues.size-1) ){
        lambdaA *= 1 / (1 + eigenA.eigenvalues[i])
    }

    val eigenB : EigenvalueDecomposition =  E.inverse().mult(hB).eigenvalueDecomposition
    var lambdaB  = 1.0
    for(i in 0..(eigenB.eigenvalues.size-1) ){
        lambdaB *= 1 / (1 + eigenB.eigenvalues[i])
    }

    val eigenD : EigenvalueDecomposition =  E.inverse().mult(hD).eigenvalueDecomposition
    var lambdaD  = 1.0
    for(i in 0..(eigenD.eigenvalues.size-1) ){
        lambdaD *= 1 / (1 + eigenD.eigenvalues[i])
    }

    val lenA = data["gender"].asStrings().distinctBy {it.hashCode()}.size.toDouble()
    val lenB = data["level"].asStrings().distinctBy {it.hashCode()}.size.toDouble()
    val lenD = data["interaction"].asStrings().distinctBy {it.hashCode()}.size.toDouble()

    var p = Y.numberOfColumns().toDouble()

    val vhA = (1 + lenA -1 + lenB - 1) - lenB
    val tA:Double
    if((p.pow(2) + vhA.pow(2) - 5) > 0){
        tA = Math.sqrt((p.pow(2)*vhA.pow(2) - 4)/(p.pow(2) + vhA.pow(2) -5))
    }else{tA = 1.0}

    val vhB= (1 + lenA -1 + lenB - 1) - lenA
    val tB:Double
    if((p.pow(2) + vhB.pow(2) - 5) > 0){
        tB = Math.sqrt((p.pow(2)*vhB.pow(2) - 4)/(p.pow(2) + vhB.pow(2) -5))
    }else{tB = 1.0}

    val vhD= (1 + lenD -1) - (lenA - 1) - lenB
    val tD:Double
    if((p.pow(2) + vhD.pow(2) - 5) > 0){
        tD = Math.sqrt((p.pow(2)*vhD.pow(2) - 4)/(p.pow(2) + vhD.pow(2) -5))
    }else{tD = 1.0}

    val ve =Y.numberOfRows().toDouble() - vhA - vhB - vhD - 1

    val rA = ve - (p-vhA+1)/2
    val rB = ve - (p-vhB+1)/2
    val rD = ve - (p-vhD+1)/2

    val fA= (p*vhA-2)/4
    val fB= (p*vhB-2)/4
    val fD= (p*vhD-2)/4

    val fCalcA = ((1-lambdaA.pow(1/tA))/(lambdaA.pow(1/tA)))*((rA*tA - 2*fA)/(p*vhA))
    val fCalcB = ((1-lambdaB.pow(1/tB))/(lambdaB.pow(1/tB)))*((rB*tB - 2*fB)/(p*vhB))
    val fCalcD = ((1-lambdaD.pow(1/tD))/(lambdaD.pow(1/tD)))*((rD*tD - 2*fD)/(p*vhD))

    val pvalueA : Double = LinearModelUtils.Ftest(fCalcA, p*vhA, (rA*tA)-(2*fA))
    val pvalueB : Double = LinearModelUtils.Ftest(fCalcB, p*vhB, (rB*tB)-(2*fB))
    val pvalueD : Double = LinearModelUtils.Ftest(fCalcD, p*vhD, (rD*tD)-(2*fD))

    println("$pvalueA, $pvalueB, $pvalueD")

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