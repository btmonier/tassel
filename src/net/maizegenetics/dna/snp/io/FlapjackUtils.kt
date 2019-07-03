@file:JvmName("FlapjackUtils")

package net.maizegenetics.dna.snp.io

import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.dna.snp.ImportUtils
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants.UNDEFINED_ALLELE_STR
import net.maizegenetics.util.LoggingUtils
import net.maizegenetics.util.Utils
import org.apache.log4j.Logger


/**
 * @author Terry Casstevens
 * Created July 02, 2019
 */

private val myLogger = Logger.getLogger("net.maizegenetics.dna.snp.io.FlapjackUtils")

private const val FLAPJACK_MISSING = "-"

// Translation from TASSEL byte encoding to Flapjack character.
// Differences from TASSEL.  N (Unknown), - (Deletion), and + (Insertion)
// are converted to FLAPJACK_MISSING (-)
// UNDEFINED_ALLELE_STR should never end up in the file.
private val FLAPJACK_CHARS = arrayOf("A", "C", "G", "T", FLAPJACK_MISSING, FLAPJACK_MISSING,
        UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR,
        UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR,
        UNDEFINED_ALLELE_STR, FLAPJACK_MISSING)


/**
 * Writes given genotype table to Flapjack format
 */
fun writeToFlapjack(genotypes: GenotypeTable, filename: String, delimChar: Char = '\t'): String {

    if (delimChar != ' ' && delimChar != '\t') {
        throw IllegalArgumentException("FlapjackUtils: writeToFlapjack: Delimiter character must be either a blank space or a tab.")
    }

    val mapFileName = Utils.addSuffixIfNeeded(filename, ".flpjk.map")
    val genoFileName = Utils.addSuffixIfNeeded(filename, ".flpjk.geno")

    try {

        Utils.getBufferedWriter(genoFileName).use { genoWriter ->

            Utils.getBufferedWriter(mapFileName).use { mapWriter ->

                for (site in 0 until genotypes.numberOfSites()) {
                    mapWriter.write(genotypes.siteName(site))
                    mapWriter.write(delimChar.toInt())
                    mapWriter.write(genotypes.chromosomeName(site))
                    mapWriter.write(delimChar.toInt())
                    mapWriter.write(Integer.toString(genotypes.chromosomalPosition(site)))
                    mapWriter.write("\n")
                    genoWriter.write(delimChar.toInt())
                    genoWriter.write(genotypes.siteName(site))
                }

            }

            genoWriter.write("\n")

            for (taxa in 0 until genotypes.numberOfTaxa()) {

                genoWriter.write(genotypes.taxaName(taxa))

                for (site in 0 until genotypes.numberOfSites()) {

                    genoWriter.write(delimChar.toInt())

                    val alleles = genotypes.genotypeArray(taxa, site)
                    when (alleles.size) {
                        1 -> {
                            genoWriter.write(FLAPJACK_CHARS[alleles[0].toInt()])
                        }
                        2 -> {
                            genoWriter.write(FLAPJACK_CHARS[alleles[0].toInt()])
                            genoWriter.write("/")
                            genoWriter.write(FLAPJACK_CHARS[alleles[1].toInt()])
                        }
                        else -> {
                            genoWriter.write(FLAPJACK_CHARS[alleles[0].toInt()])
                            genoWriter.write("/")
                            genoWriter.write(FLAPJACK_CHARS[alleles[1].toInt()])
                            genoWriter.write("/")
                            genoWriter.write(FLAPJACK_CHARS[alleles[2].toInt()])
                        }
                    }

                }

                genoWriter.write("\n")
            }
        }

        return "$mapFileName and $genoFileName"

    } catch (e: Exception) {
        myLogger.debug(e.message, e)
        throw IllegalStateException("Problem writing Flapjack files: $mapFileName and $genoFileName\n${e.message}")
    }

}

fun main() {
    LoggingUtils.setupDebugLogging()
    val filename = "/Users/terry/terry/tassel-5-source/mdp_genotype.hmp.txt"
    val genotype = ImportUtils.read(filename)
    writeToFlapjack(genotype, "/Users/terry/terry/tassel-5-source/mdp_genotype")
}