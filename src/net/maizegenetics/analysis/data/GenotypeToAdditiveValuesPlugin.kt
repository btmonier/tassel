package net.maizegenetics.analysis.data

import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.dna.snp.GenotypeTableUtils
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.util.ColumnMatrix
import java.awt.Frame
import javax.swing.ImageIcon

/**
 * @author Terry Casstevens
 * Created December 06, 2018
 */

class GenotypeToAdditiveValuesPlugin(parentFrame: Frame?, isInteractive: Boolean) : AbstractPlugin(parentFrame, isInteractive) {

    override fun processData(input: DataSet?): DataSet {

        val temp = input?.getDataOfType(GenotypeTable::class.java)
        if (temp?.size != 1) {
            throw IllegalArgumentException("GenotypeToAdditiveValuesPlugin: processData: must input a genotype")
        }
        val genotype = temp[0].data as GenotypeTable

        val result = ColumnMatrix.Builder(genotype.numberOfTaxa(), genotype.numberOfSites())

        for (s in 0 until genotype.numberOfSites()) {
            val site = genotype.genotypeAllTaxa(s)
            val alleleCounts = AlleleFreqCache.allelesSortedByFrequencyNucleotide(site)
            val majorAllele = AlleleFreqCache.majorAllele(alleleCounts)
            for (t in 0 until genotype.numberOfTaxa()) {
                var value: Byte = 0
                val alleles = GenotypeTableUtils.getDiploidValues(site[t])
                if (alleles[0] != majorAllele) value++
                if (alleles[1] != majorAllele) value++
                result.set(t, s, value)
            }
        }

        return DataSet(Datum("${temp[0].name} Additive Values", result.build(), null), this)

    }

    fun runPlugin(genotype: GenotypeTable): ColumnMatrix {
        return performFunction(DataSet.getDataSet(genotype)).getDataOfType(ColumnMatrix::class.java)[0].data as ColumnMatrix
    }

    override fun getIcon(): ImageIcon? {
        return null
    }

    override fun getButtonName(): String {
        return "Convert Genotype to Additive Values"
    }

    override fun getToolTipText(): String {
        return "Convert Genotype to Additive Values"
    }

}
