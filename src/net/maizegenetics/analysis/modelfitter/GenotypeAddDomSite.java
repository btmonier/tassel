package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.stats.linearmodels.ModelEffect;

public class GenotypeAddDomSite extends GenotypeAdditiveSite implements DominantSite  {

    public GenotypeAddDomSite(int site, String chr, int pos, String id,
                              CRITERION selectionCriterion, byte[] genotype,
                              byte majorAllele, double majorAlleleFrequency) {

        super(site, chr,pos, id, selectionCriterion, genotype, majorAllele, majorAlleleFrequency);
    }

    @Override
    public ModelEffect getModelEffect() {
        if (taxaIndex == null)
            return getModelEffectNoReindex();
        else
            return getModelEffectWithReindex();
    }

    @Override
    public ModelEffect getModelEffect(int[] subset) {
        if (taxaIndex == null)
            return getModelEffectNoReindex(subset);
        else
            return getModelEffectWithReindex(subset);
    }

    private ModelEffect getModelEffectNoReindex() {
        //Use ModelEffectUtils method to generate integer levels
        //use integer levels to instantiate a Factor(?)ModelEffect
        //the question is how to deal with missing data
        return null;
    }

    private ModelEffect getModelEffectWithReindex() {
        return null;
    }

    private ModelEffect getModelEffectNoReindex(int[] subset) {
        return null;
    }

    private ModelEffect getModelEffectWithReindex(int[] subset) {
        return null;
    }
}
