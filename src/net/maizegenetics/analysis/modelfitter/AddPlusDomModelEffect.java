package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import org.apache.log4j.Logger;

import java.util.Arrays;

public class AddPlusDomModelEffect implements ModelEffect {

    private Object id;
    private double[] addCovariate;
    private double[] domCovariate;
    private CovariateModelEffect addModelEffect;
    private CovariateModelEffect domModelEffect;
    private double minHets = 5;

    public AddPlusDomModelEffect(Object id, AdditiveSite addSite) {
        this(id, addSite.getCovariate());
    }

    public AddPlusDomModelEffect(Object id, double[] additiveCovariate) {
        addCovariate = additiveCovariate;
        domCovariate = Arrays.stream(addCovariate).map(add -> 1.0 - Math.abs(add - 1)).toArray();
        double domSum = Arrays.stream(domCovariate).sum();

        addModelEffect = new CovariateModelEffect(addCovariate, id);
        if (domSum >= minHets) {
            domModelEffect = new CovariateModelEffect(domCovariate, id);
        } else {
            domModelEffect = null;
        }
    }

    @Override
    public Object getID() {
        return id;
    }

    @Override
    public void setID(Object id) {
        this.id = id;
    }

    @Override
    public int getSize() {
        return addCovariate.length;
    }

    @Override
    public DoubleMatrix getX() {
        if (domModelEffect == null) {
            return addModelEffect.getX();
        } else {
            return addModelEffect.getX().concatenate(domModelEffect.getX(), false);
        }

    }

    @Override
    public DoubleMatrix getXtX() {
        return getX().crossproduct();
    }

    @Override
    public DoubleMatrix getXty(double[] y) {
        return getX().crossproduct(DoubleMatrixFactory.DEFAULT.make(1, y.length, y));
    }

    @Override
    public DoubleMatrix getyhat(DoubleMatrix beta) {
        return getX().crossproduct(beta);
    }

    @Override
    public DoubleMatrix getyhat(double[] beta) {
        return getyhat(DoubleMatrixFactory.DEFAULT.make(1, beta.length, beta));
    }

    @Override
    public int[] getLevelCounts() {
        if (domModelEffect == null) {
            return new int[]{addCovariate.length};
        } else return new int[]{addCovariate.length, domCovariate.length};
    }

    @Override
    public int getNumberOfLevels() {
        if (domModelEffect == null) {
            return 1;
        } else return 2;

    }

    @Override
    public int getEffectSize() {
        return getNumberOfLevels();
    }

    @Override
    public ModelEffect getCopy() {
        return new AddPlusDomModelEffect(id, addCovariate);
    }

    @Override
    public ModelEffect getSubSample(int[] sample) {
        double[] subsample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
            subsample[i] = addCovariate[sample[i]];
        }
        return new AddPlusDomModelEffect(id, subsample);
    }
}
