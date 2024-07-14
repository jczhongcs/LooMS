package edu.uw.waterlooms.peptideMatch;

import org.apache.commons.math3.util.FastMath;

public class SimilarityCalculate {

    /** sum of x values */
    private double sumX = 0d;

    /** total variation in x (sum of squared deviations from xbar) */
    private double sumXX = 0d;

    /** sum of y values */
    private double sumY = 0d;

    /** total variation in y (sum of squared deviations from ybar) */
    private double sumYY = 0d;

    /** sum of products */
    private double sumXY = 0d;

    /** number of observations */
    private long n = 0;

    /** mean of accumulated x values, used in updating formulas */
    private double xbar = 0;

    /** mean of accumulated y values, used in updating formulas */
    private double ybar = 0;

    /** include an intercept or not */
    private final boolean hasIntercept;

    public SimilarityCalculate(boolean hasIntercept) {
        this.hasIntercept = hasIntercept;
    }
    public void clear() {
        sumX = 0d;
        sumXX = 0d;
        sumY = 0d;
        sumYY = 0d;
        sumXY = 0d;
        n = 0;
    }
    public void addData(final double x, final double y) {

        if (n == 0) {
            xbar = x;
            ybar = y;
        } else {
            if( hasIntercept ){
                final double fact1 = 1.0 + n;
                final double fact2 = n / (1.0 + n);
                final double dx = x - xbar;
                final double dy = y - ybar;
                sumXX += dx * dx * fact2;
                sumYY += dy * dy * fact2;
                sumXY += dx * dy * fact2;
                xbar += dx / fact1;
                ybar += dy / fact1;
            }
        }
        if( !hasIntercept ){
            sumXX += x * x ;
            sumYY += y * y ;
            sumXY += x * y ;
        }
        sumX += x;
        sumY += y;
        n++;
        return;
    }
    public double getCosSimilarity()
    {
        if (FastMath.abs(sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(sumYY) < 10 * Double.MIN_VALUE) {
            return 0.0;
        }else
        {
            return sumXY/FastMath.sqrt(sumXX*sumYY);
        }
    }
    public double getSlope() {
        if (n < 2) {
            return Double.NaN; //not enough data
        }
        if (FastMath.abs(sumXX) < 10 * Double.MIN_VALUE) {
            return Double.NaN; //not enough variation in x
        }
        return sumXY / sumXX;
    }
    public double getRSquare() {
        double ssto = getTotalSumSquares();
        return (ssto - getSumSquaredErrors()) / ssto;
    }
    public double getTotalSumSquares() {
        if (n < 2) {
            return Double.NaN;
        }
        return sumYY;
    }
    public double getSumSquaredErrors() {
        return FastMath.max(0d, sumYY - sumXY * sumXY / sumXX);
    }
    public double getR() {
        double b1 = getSlope();
        double result = FastMath.sqrt(getRSquare());
        if (b1 < 0) {
            result = -result;
        }
        return result;
    }
}
