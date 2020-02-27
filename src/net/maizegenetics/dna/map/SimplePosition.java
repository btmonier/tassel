package net.maizegenetics.dna.map;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.util.GeneralAnnotation;

/**
 * @author Terry Casstevens Created August 29, 2019
 */
public class SimplePosition extends GeneralPosition {

    public SimplePosition(Chromosome chromosome, int position) {
        super(chromosome, position);
    }

    @Override
    public String[] getKnownVariants() {
        return new String[0];
    }

    @Override
    public float getGlobalMAF() {
        return Float.NaN;
    }

    @Override
    public float getGlobalSiteCoverage() {
        return Float.NaN;
    }

    @Override
    public byte getAllele(WHICH_ALLELE alleleType) {
        return GenotypeTable.UNKNOWN_ALLELE;
    }

    @Override
    public GeneralAnnotation getAnnotation() {
        return null;
    }

    @Override
    public int compareTo(Position o) {
        if (o instanceof SimplePosition) {
            int result = myChromosome.compareTo(o.getChromosome());
            if (result != 0) {
                return result;
            }
            return Integer.compare(myPosition, o.getPosition());
        } else {
            return super.compareTo(o);
        }
    }
}
