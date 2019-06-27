package net.maizegenetics.analysis.association

import com.google.common.collect.HashMultiset
import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition
import net.maizegenetics.phenotype.GenotypePhenotype
import net.maizegenetics.phenotype.Phenotype
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.PluginParameter
import net.maizegenetics.stats.linearmodels.FactorModelEffect
import net.maizegenetics.stats.linearmodels.LinearModelUtils
import net.maizegenetics.stats.linearmodels.ModelEffect
import net.maizegenetics.stats.linearmodels.ModelEffectUtils
import net.maizegenetics.util.TableReportBuilder
import org.apache.log4j.Logger
import java.awt.Frame
import java.util.*
import javax.swing.ImageIcon
import kotlin.math.pow
import kotlin.random.Random

/**
 * @author Samuel Fernandes and Terry Casstevens
 * Created June 24, 2019
 */

class ManovaPlugin(parentFrame: Frame?, isInteractive: Boolean) : AbstractPlugin(parentFrame, isInteractive) {

    private val myLogger = Logger.getLogger(ManovaPlugin::class.java)

    private var usePermutations = PluginParameter.Builder("usePerm", true, Boolean::class.javaObjectType)
            .description("Should permutations be used to set the enter and exit limits for stepwise regression? A permutation test will be used to determine the enter limit. The exit limit will be set to 2 times the enter limit.")
            .guiName("Use permutations")
            .build()

    private var numberOfPermutations = PluginParameter.Builder("nPerm", 1000, Int::class.javaObjectType)
            .description("The number of permutations used to determine the enter limit.")
            .guiName("Number of permutations")
            .dependentOnParameter(usePermutations)
            .build()

    private var enterLimit = PluginParameter.Builder("enterLimit", 1e-5, Double::class.javaObjectType)
            .description("When p-value is the model selection criteria, model fitting will stop when the next term chosen has a p-value greater than the enterLimit. This value will be over-ridden the permutation test, if used.")
            .guiName("enterLimit")
            .dependentOnParameter(usePermutations, false)
            .build()

    private var exitLimit = PluginParameter.Builder("exitLimit", 2e-5, Double::class.javaObjectType)
            .description("During the backward step of model fitting if p-value has been chosen as the selection criterion, if the term in model with the highest p-value has a p-value > exitLimit, it will be removed from the model.")
            .guiName("exitLimit")
            .dependentOnParameter(usePermutations, false)
            .build()

    private var isNested = PluginParameter.Builder("isNested", true, Boolean::class.javaObjectType)
            .description("Should SNPs/markers be nested within a factor, such as family?")
            .guiName("")
            .build()

    private var nestingFactor = PluginParameter.Builder("nestFactor", null, String::class.java)
            .guiName("Nesting factor")
            .description("Nest markers within this factor. This parameter cannot be set from the command line. Instead, the first factor in the data set will be used.")
            .dependentOnParameter(isNested)
            .objectListSingleSelect()
            .build()

    private var myGenotypeTable: PluginParameter<GenotypeTable.GENOTYPE_TABLE_COMPONENT> = PluginParameter.Builder("genotypeComponent", GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT::class.java)
            .genotypeTable()
            .range(GenotypeTable.GENOTYPE_TABLE_COMPONENT.values())
            .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
            .build()

    private var createAnova = PluginParameter.Builder("anova", true, Boolean::class.javaObjectType)
            .description("Create pre- and post-scan anova reports.")
            .guiName("Create anova reports")
            .build()

    private var createEffects = PluginParameter.Builder("effects", true, Boolean::class.javaObjectType)
            .description("Create a report of marker effects based on the scan results.")
            .guiName("Create effects report")
            .build()

    private var createStep = PluginParameter.Builder("step", true, Boolean::class.javaObjectType)
            .description("Create a report of the which markers enter and leave the model as it is being fit.")
            .guiName("Create step report")
            .build()

    private var createResiduals = PluginParameter.Builder("residuals", false, Boolean::class.javaObjectType)
            .description("Create a phenotype dataset of model residuals for each chromosome. For each chromosome, the residuals will be calculated from a model with all terms EXCEPT the markers on that chromosome.")
            .guiName("Create residuals")
            .build()

    private var writeFiles = PluginParameter.Builder("saveToFile", false, Boolean::class.javaObjectType)
            .description("Should the requested output be written to files?")
            .guiName("Write to files")
            .build()

    private var outputName = PluginParameter.Builder("savePath", "", String::class.java)
            .description("The base file path for the save files. Each file saved will add a descriptive name to the base name.")
            .guiName("Base file path")
            .outFile()
            .dependentOnParameter(writeFiles)
            .build()

    private var maximumNumberOfVariantsInModel = PluginParameter.Builder("maxQTN", 100, Int::class.javaObjectType)
            .description("maximum number of QTN to be fit in the model")
            .guiName("Maximum QTN Number")
            .build()

    private lateinit var myGenoPheno: GenotypePhenotype
    private lateinit var myDatasetName: String
    private val myFactorNameList: MutableList<String> = ArrayList()

    //TableReport builders
    private val manovaReportBuilder =
            TableReportBuilder.getInstance("Manova", arrayOf("SiteID", "Chr", "Position", "action", "approx_F", "num_df", "den_df", "probF"))
    private val permutationReportBuilder =
            TableReportBuilder.getInstance("Empirical Null", arrayOf("Trait", "p-value"))
    private val stepsReportBuilder =
            TableReportBuilder.getInstance("Steps", arrayOf("SiteID", "Chr", "Position", "action", "approx_F", "num_df", "den_df", "probF"))

    override fun preProcessParameters(input: DataSet?) {

        DoubleMatrixFactory.setDefault(DoubleMatrixFactory.FactoryType.ejml)

        //input data should be a single GenotypePhenotype
        val datumList = input?.getDataOfType(GenotypePhenotype::class.java)
        if (datumList!!.size != 1)
            throw IllegalArgumentException("Choose exactly one dataset that has combined genotype and phenotype data.")
        myGenoPheno = datumList[0].data as GenotypePhenotype
        myDatasetName = datumList[0].name

        myGenoPheno.phenotype().attributeListOfType(Phenotype.ATTRIBUTE_TYPE.factor).stream()
                .map { pa -> pa.name() }
                .forEach { myFactorNameList.add(it) }

        myFactorNameList.add("None")

        if (myFactorNameList.isEmpty()) myFactorNameList.add("None")
        nestingFactor = PluginParameter<String>(nestingFactor, myFactorNameList)

    }

    override fun processData(input: DataSet): DataSet? {

        var xR = DoubleMatrixFactory.DEFAULT.make(myGenoPheno.numberOfObservations(), 1, 1.0)
        val Y = createY()

        val modelEffectList = ArrayList<ModelEffect>()
        val step1 = forwardStep(Y, xR, modelEffectList, null)
        xR = step1.first
        var snpInModel = step1.second
        while (xR != null && modelEffectList.size <= maximumNumberOfVariantsInModel()) {
            var step2 = forwardStep(Y, xR, modelEffectList, snpInModel)
            xR = step2.first
            snpInModel = step2.second
        }

        return null
    }

    fun forwardStep(Y: DoubleMatrix, xR: DoubleMatrix, modelEffectList: MutableList<ModelEffect>, snpInModel: MutableList<String>?): Pair<DoubleMatrix?, MutableList<String>?> {
        val nSites = myGenoPheno.genotypeTable().numberOfSites()
        var minPval = 1.0
        var bestModelEffect: ModelEffect? = null

        //TODO do not test sites already in modelEffectList
        for (sitenum in 0 until nSites) {
            if (snpInModel == null || !snpInModel.contains(myGenoPheno.genotypeTable().siteName(sitenum))) {
                val genotypesForSite = imputeNsInGenotype(myGenoPheno.getStringGenotype(sitenum))
                val modelEffect = FactorModelEffect(ModelEffectUtils.getIntegerLevels(genotypesForSite),
                        true, myGenoPheno.genotypeTable().siteName(sitenum))
                val siteDesignMatrix = xR.concatenate(modelEffect.x, false)
                println(xR)
                val pval = manovaPvalue(Y, siteDesignMatrix, xR)
                if (pval < minPval) {
                    minPval = pval
                    bestModelEffect = modelEffect
                }
            }
        }

        if (minPval <= enterLimit() && bestModelEffect != null) {
            modelEffectList.add(bestModelEffect)
            println("${bestModelEffect.id}, pval = $minPval")
            snpInModel?.add(bestModelEffect.id.toString())
            return Pair(xR.concatenate(bestModelEffect.x, false), snpInModel)
        }
        return Pair(null, null)
    }

    fun imputeNsInGenotype(genotypes: Array<String>): Array<String> {
        //TODO write test to make sure this works
        val alleleCounter = HashMultiset.create<String>()
        for (allele in genotypes) if (!allele.equals("N")) alleleCounter.add(allele)
        val totalCount = alleleCounter.count().toDouble()
        val alleleProportionList = ArrayList<Pair<Double, String>>()
        var cumulativeSum = 0
        for (entry in alleleCounter.entrySet()) {
            cumulativeSum += entry.count
            alleleProportionList.add(Pair(cumulativeSum / totalCount, entry.element))
        }

        return genotypes.map {
            if (it.equals("N")) {
                val ran = Random.nextDouble()
                var ndx = 0
                while (ran >= alleleProportionList[ndx].first) ndx++
                alleleProportionList[ndx].second
            } else it
        }.toTypedArray()

    }

    fun createY(): DoubleMatrix {
        val dataAttributeList = myGenoPheno.phenotype().attributeListOfType(Phenotype.ATTRIBUTE_TYPE.data)
        val nTraits = dataAttributeList.size
        val nObs = myGenoPheno.phenotype().numberOfObservations()
        val Y = DoubleMatrixFactory.DEFAULT.make(nObs, nTraits)
        for (t in 0 until nTraits) {
            val traitvals = dataAttributeList[t].allValues() as FloatArray
            traitvals.forEachIndexed { index, fl ->
                val traitvalue = fl.toDouble()
                if (traitvalue.isNaN()) throw java.lang.IllegalArgumentException("A trait value is missing. No missing trait values are allowed.")
                Y.set(index, t, traitvalue)
            }
        }
        return Y
    }


    fun calcBeta(Y: DoubleMatrix, X: DoubleMatrix): BetaValue {
        val XtX: DoubleMatrix = X.crossproduct()
        val XtY: DoubleMatrix = X.crossproduct(Y)
        val XtXinv: DoubleMatrix = XtX.generalizedInverse()
        val B: DoubleMatrix = XtXinv.mult(XtY)
        val H: DoubleMatrix = B.transpose().mult(XtY)
        return BetaValue(B, H)
    }

    /**
     * Y = Phenotypic data
     * xF = Full model
     * xR = Reduced model
     * returns Wilk's lambda pvalue
     */
    private fun manovaPvalue(Y: DoubleMatrix, xF: DoubleMatrix, xR: DoubleMatrix): Double {
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
        if (df == 0.0) {
            return 1.0
        }
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

    /**
     * Y = Phenotypic data
     * xI = intercept
     * returns Wilk's lambda pvalue
     */
    private fun manovaPvalue(Y: DoubleMatrix, xI: DoubleMatrix): Double {
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

    override fun getIcon(): ImageIcon? {
        return null
    }

    override fun getButtonName(): String {
        return "Multi-trait Stepwise"
    }

    override fun getToolTipText(): String {
        return "Multi-trait Stepwise"
    }

    override fun getCitation(): String {
        return "reference..."
    }

    override fun pluginUserManualURL(): String {
        return "https://bitbucket.org/tasseladmin/tassel­5­source/wiki/UserManual/..."
    }


    /**
     * Should permutations be used to set the enter and exit
     * limits for stepwise regression? A permutation test
     * will be used to determine the enter limit. The exit
     * limit will be set to 2 times the enter limit.
     *
     * @return Use permutations
     */
    fun usePermutations(): Boolean {
        return usePermutations.value()
    }

    /**
     * Set Use permutations. Should permutations be used to
     * set the enter and exit limits for stepwise regression?
     * A permutation test will be used to determine the enter
     * limit. The exit limit will be set to 2 times the enter
     * limit.
     *
     * @param value Use permutations
     *
     * @return this plugin
     */
    fun usePermutations(value: Boolean): ManovaPlugin {
        usePermutations = PluginParameter(usePermutations, value)
        return this
    }

    /**
     * The number of permutations used to determine the enter
     * limit.
     *
     * @return Number of permutations
     */
    fun numberOfPermutations(): Int {
        return numberOfPermutations.value()
    }

    /**
     * Set Number of permutations. The number of permutations
     * used to determine the enter limit.
     *
     * @param value Number of permutations
     *
     * @return this plugin
     */
    fun numberOfPermutations(value: Int): ManovaPlugin {
        numberOfPermutations = PluginParameter(numberOfPermutations, value)
        return this
    }

    /**
     * When p-value is the model selection criteria, model
     * fitting will stop when the next term chosen has a p-value
     * greater than the enterLimit. This value will be over-ridden
     * the permutation test, if used.
     *
     * @return enterLimit
     */
    fun enterLimit(): Double {
        return enterLimit.value()
    }

    /**
     * Set enterLimit. When p-value is the model selection
     * criteria, model fitting will stop when the next term
     * chosen has a p-value greater than the enterLimit. This
     * value will be over-ridden the permutation test, if
     * used.
     *
     * @param value enterLimit
     *
     * @return this plugin
     */
    fun enterLimit(value: Double): ManovaPlugin {
        enterLimit = PluginParameter(enterLimit, value)
        return this
    }

    /**
     * During the backward step of model fitting if p-value
     * has been chosen as the selection criterion, if the
     * term in model with the highest p-value has a p-value
     * > exitLimit, it will be removed from the model.
     *
     * @return exitLimit
     */
    fun exitLimit(): Double {
        return exitLimit.value()
    }

    /**
     * Set exitLimit. During the backward step of model fitting
     * if p-value has been chosen as the selection criterion,
     * if the term in model with the highest p-value has a
     * p-value > exitLimit, it will be removed from the model.
     *
     * @param value exitLimit
     *
     * @return this plugin
     */
    fun exitLimit(value: Double): ManovaPlugin {
        exitLimit = PluginParameter(exitLimit, value)
        return this
    }

    /**
     * Should SNPs/markers be nested within a factor, such
     * as family?
     *
     * @return Is Nested
     */
    fun isNested(): Boolean {
        return isNested.value()
    }

    /**
     * Set Is Nested. Should SNPs/markers be nested within
     * a factor, such as family?
     *
     * @param value Is Nested
     *
     * @return this plugin
     */
    fun isNested(value: Boolean): ManovaPlugin {
        isNested = PluginParameter(isNested, value)
        return this
    }

    /**
     * Nest markers within this factor. This parameter cannot
     * be set from the command line. Instead, the first factor
     * in the data set will be used.
     *
     * @return Nesting factor
     */
    fun nestingFactor(): String {
        return nestingFactor.value()
    }

    /**
     * Set Nesting factor. Nest markers within this factor.
     * This parameter cannot be set from the command line.
     * Instead, the first factor in the data set will be used.
     *
     * @param value Nesting factor
     *
     * @return this plugin
     */
    fun nestingFactor(value: String): ManovaPlugin {
        nestingFactor = PluginParameter(nestingFactor, value)
        return this
    }

    /**
     * If the genotype table contains more than one type of
     * genotype data, choose the type to use for the analysis.
     *
     * @return Genotype Component
     */
    fun genotypeTable(): GenotypeTable.GENOTYPE_TABLE_COMPONENT {
        return myGenotypeTable.value()
    }

    /**
     * Set Genotype Component. If the genotype table contains
     * more than one type of genotype data, choose the type
     * to use for the analysis.
     *
     * @param value Genotype Component
     *
     * @return this plugin
     */
    fun genotypeTable(value: GenotypeTable.GENOTYPE_TABLE_COMPONENT): ManovaPlugin {
        myGenotypeTable = PluginParameter(myGenotypeTable, value)
        return this
    }

    /**
     * Create pre- and post-scan anova reports.
     *
     * @return Create anova reports
     */
    fun createAnova(): Boolean {
        return createAnova.value()
    }

    /**
     * Set Create anova reports. Create pre- and post-scan
     * anova reports.
     *
     * @param value Create anova reports
     *
     * @return this plugin
     */
    fun createAnova(value: Boolean): ManovaPlugin {
        createAnova = PluginParameter(createAnova, value)
        return this
    }

    /**
     * Create a report of marker effects based on the scan
     * results.
     *
     * @return Create effects report
     */
    fun createEffects(): Boolean {
        return createEffects.value()
    }

    /**
     * Set Create effects report. Create a report of marker
     * effects based on the scan results.
     *
     * @param value Create effects report
     *
     * @return this plugin
     */
    fun createEffects(value: Boolean): ManovaPlugin {
        createEffects = PluginParameter(createEffects, value)
        return this
    }

    /**
     * Create a report of the which markers enter and leave
     * the model as it is being fit.
     *
     * @return Create step report
     */
    fun createStep(): Boolean {
        return createStep.value()
    }

    /**
     * Set Create step report. Create a report of the which
     * markers enter and leave the model as it is being fit.
     *
     * @param value Create step report
     *
     * @return this plugin
     */
    fun createStep(value: Boolean): ManovaPlugin {
        createStep = PluginParameter(createStep, value)
        return this
    }

    /**
     * Create a phenotype dataset of model residuals for each
     * chromosome. For each chromosome, the residuals will
     * be calculated from a model with all terms EXCEPT the
     * markers on that chromosome.
     *
     * @return Create residuals
     */
    fun createResiduals(): Boolean {
        return createResiduals.value()
    }

    /**
     * Set Create residuals. Create a phenotype dataset of
     * model residuals for each chromosome. For each chromosome,
     * the residuals will be calculated from a model with
     * all terms EXCEPT the markers on that chromosome.
     *
     * @param value Create residuals
     *
     * @return this plugin
     */
    fun createResiduals(value: Boolean): ManovaPlugin {
        createResiduals = PluginParameter(createResiduals, value)
        return this
    }

    /**
     * Should the requested output be written to files?
     *
     * @return Write to files
     */
    fun writeFiles(): Boolean {
        return writeFiles.value()
    }

    /**
     * Set Write to files. Should the requested output be
     * written to files?
     *
     * @param value Write to files
     *
     * @return this plugin
     */
    fun writeFiles(value: Boolean): ManovaPlugin {
        writeFiles = PluginParameter(writeFiles, value)
        return this
    }

    /**
     * The base file path for the save files. Each file saved
     * will add a descriptive name to the base name.
     *
     * @return Base file path
     */
    fun outputName(): String {
        return outputName.value()
    }

    /**
     * Set Base file path. The base file path for the save
     * files. Each file saved will add a descriptive name
     * to the base name.
     *
     * @param value Base file path
     *
     * @return this plugin
     */
    fun outputName(value: String): ManovaPlugin {
        outputName = PluginParameter(outputName, value)
        return this
    }

    /**
     * maximum number of QTN to be fit in the model
     *
     * @return Maximum QTN Number
     */
    fun maximumNumberOfVariantsInModel(): Int {
        return maximumNumberOfVariantsInModel.value()
    }

    /**
     * Set Maximum QTN Number. maximum number of QTN to be
     * fit in the model
     *
     * @param value Maximum QTN Number
     *
     * @return this plugin
     */
    fun maximumNumberOfVariantsInModel(value: Int): ManovaPlugin {
        maximumNumberOfVariantsInModel = PluginParameter(maximumNumberOfVariantsInModel, value)
        return this
    }

}

//fun main(args : Array<String>) {
//    GeneratePluginCode.generate(ManovaPlugin::class.java)
//}

data class BetaValue(val B: DoubleMatrix, val H: DoubleMatrix)