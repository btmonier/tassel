package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.stats.linearmodels.ModelEffect;

public interface DominantSite {
    /**
     * @return		the ModelEffect for this site
     */
    ModelEffect getModelEffect();


    /**
     * @param subset	an int array indexing a subset of taxa
     * @return			the ModelEffect for the subset of taxa for this sites
     */
    ModelEffect getModelEffect(int[] subset);
}
