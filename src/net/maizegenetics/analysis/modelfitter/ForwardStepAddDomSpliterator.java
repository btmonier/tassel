package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.stats.linearmodels.ModelEffect;

import java.util.List;
import java.util.function.Consumer;

public class ForwardStepAddDomSpliterator extends ForwardStepAdditiveSpliterator {

    public ForwardStepAddDomSpliterator(List<AdditiveSite> siteList, List<ModelEffect> baseModel,
                                        double[] y) {
        super(siteList, baseModel, y);
    }

    //implements an additive plus dominance model
    @Override
    public boolean tryAdvance(Consumer<? super AdditiveSite> action) {
        if (origin == end)
            return false;
        AdditiveSite as = mySites.get(origin);

        switch (as.selectionCriterion()) {
            case pval:
                AddPlusDomModelEffect me = new AddPlusDomModelEffect(as, as.getCovariate());
                plm.testNewModelEffect(me);
                as.criterionValue(-Math.log10(plm.getp())); //use -log10 p as criterion so that objective is to maximize it
                break;
            case aic:
                plm.testNewModelEffect(new AddPlusDomModelEffect(as, as.getCovariate()));
                double rss = plm.getErrorSS();
                as.criterionValue(nobs * Math.log(rss / nobs) + 2 * (baseModeldf + 1));
                break;
            case bic:
                plm.testNewModelEffect(new AddPlusDomModelEffect(as, as.getCovariate()));
                rss = plm.getErrorSS();
                as.criterionValue(nobs * Math.log(rss / nobs) + Math.log(nobs) * (baseModeldf + 1));
                break;
            case mbic:
                plm.testNewModelEffect(new AddPlusDomModelEffect(as, as.getCovariate()));
                rss = plm.getErrorSS();
                as.criterionValue(nobs * Math.log(rss / nobs) + Math.log(nobs) * (baseModeldf + 1) + 2
                        * (baseModeldf + 1) * Math.log(nsites / 2.2 - 1));
                break;
        }

        action.accept(as);
        origin++;
        return true;
    }

}
