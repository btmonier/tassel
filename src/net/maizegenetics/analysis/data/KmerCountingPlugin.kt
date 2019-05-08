package net.maizegenetics.analysis.data

import htsjdk.samtools.fastq.FastqReader
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.PluginParameter
import org.apache.log4j.Logger
import java.awt.Frame
import java.io.File
import javax.swing.ImageIcon

/**
 * @author Terry Casstevens
 * Created May 06, 2019
 */

private const val NUM_BASES_PER_HASH = 14

class KmerCountingPlugin(parentFrame: Frame?, isInteractive: Boolean) : AbstractPlugin(parentFrame, isInteractive) {

    private val myLogger = Logger.getLogger(KmerCountingPlugin::class.java)

    private var fastqFile = PluginParameter.Builder("fastqFile", null, String::class.java)
            .description("Fastq filename")
            .required(true)
            .inFile()
            .build()

    private var kmerLength = PluginParameter.Builder("kmerLength", 3, Int::class.javaObjectType)
            .description("Length of Kmer")
            .build()

    private fun numHashes() = { Math.ceil(kmerLength() / NUM_BASES_PER_HASH.toDouble()).toInt() }.invoke()

    override fun processData(input: DataSet?): DataSet? {

        val countFunction: (String) -> Unit = if (numHashes() == 1) (::countSingle) else (::countTree)

        FastqReader(File(fastqFile())).use { reader ->

            while (reader.hasNext()) {

                val record = reader.next()
                val sequence = record.readString
                sequence.windowed(size = kmerLength(), step = kmerLength(), partialWindows = false)
                        .forEach { countFunction(it) }

            }

        }

        return null

    }

    private val countMap by lazy { HashMap<String, Int>() }
    private fun countSingle(kmer: String) {
        val count = countMap[kmer]
        if (count == null) countMap.put(kmer, 1) else countMap.put(kmer, count + 1)
    }

    private val countTreeMap by lazy { HashMap<String, Any>() }
    private fun countTree(kmer: String) {

        var currentMap: HashMap<String, Any> = countTreeMap
        for (i in 0 until numHashes() - 2) {
            val seq = kmer.substring(i * NUM_BASES_PER_HASH, i * NUM_BASES_PER_HASH + NUM_BASES_PER_HASH)
            var temp = currentMap[seq] as HashMap<String, Any>?
            if (temp == null) {
                temp = HashMap()
                currentMap[seq] = temp
            }
            currentMap = temp
        }

        val seqNextToLast = kmer.substring((numHashes() - 2) * NUM_BASES_PER_HASH, (numHashes() - 2) * NUM_BASES_PER_HASH + NUM_BASES_PER_HASH)
        var counts = currentMap[seqNextToLast] as HashMap<String, Int>?
        if (counts == null) {
            counts = HashMap()
            currentMap[seqNextToLast] = counts
        }

        val seqLast = kmer.substring((numHashes() - 1) * NUM_BASES_PER_HASH, Math.min((numHashes() - 1) * NUM_BASES_PER_HASH + NUM_BASES_PER_HASH, kmer.length))
        val count = currentMap[seqLast]
        if (count == null) currentMap.put(kmer, 1) else currentMap.put(kmer, count as Int + 1)

    }

    private fun seqToInts(seq: String, start: Int): Array<String>? {

        return Array(numHashes()) { i ->
            val startIndex = start + NUM_BASES_PER_HASH * i
            seq.substring(startIndex, Math.min(startIndex + NUM_BASES_PER_HASH, start + kmerLength()))
        }

    }

    fun fastqFile(): String {
        return fastqFile.value()
    }

    fun fastqFile(value: String): KmerCountingPlugin {
        fastqFile = PluginParameter(fastqFile, value)
        return this
    }

    fun kmerLength(): Int {
        return kmerLength.value()
    }

    fun knerLength(value: Int): KmerCountingPlugin {
        kmerLength = PluginParameter(kmerLength, value)
        return this
    }

    override fun getToolTipText(): String {
        return "Kmer Counting"
    }

    override fun getIcon(): ImageIcon? {
        return null
    }

    override fun getButtonName(): String {
        return "Kmer Counting"
    }

}