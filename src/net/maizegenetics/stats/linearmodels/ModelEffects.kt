package net.maizegenetics.stats.linearmodels

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix

class NestedFactorModelEffect(nestedEffect : FactorModelEffect, withinEffect : FactorModelEffect, id : Any) : ModelEffect {
    override fun getID(): Any {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun setID(id: Any?) {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getSize(): Int {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getX(): DoubleMatrix {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getXtX(): DoubleMatrix {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getXty(y: DoubleArray?): DoubleMatrix {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getyhat(beta: DoubleMatrix?): DoubleMatrix {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getyhat(beta: DoubleArray?): DoubleMatrix {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getLevelCounts(): IntArray {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getNumberOfLevels(): Int {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getEffectSize(): Int {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getCopy(): ModelEffect {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getSubSample(sample: IntArray?): ModelEffect {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }
}