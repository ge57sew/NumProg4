package ode;

/**
 * Der klassische Runge-Kutta der Ordnung 4
 *
 * @author braeckle
 *
 */
public class RungeKutta4 implements Einschrittverfahren {

    @Override
    /**
     * {@inheritDoc}
     * Bei der Umsetzung koennen die Methoden addVectors und multScalar benutzt werden.
     */
    public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
        // TODO: done

        double[] f1 = ode.auswerten(t, y_k);
        double[] k1 = multScalar(f1, delta_t);

        double[] r1 = addVectors(y_k, multScalar(k1, 0.5));
        double[] f2 = ode.auswerten(t+delta_t/2.0, r1);
        double[] k2 = multScalar(f2, delta_t);

        double[] r2 = addVectors(y_k, multScalar(k2, 0.5));
        double[] f3 = ode.auswerten(t+delta_t/2.0, r2);
        double[] k3 = multScalar(f3, delta_t);

        double[] r3 = addVectors(y_k, k3);
        double[] f4 = ode.auswerten(t+delta_t  , r3);
        double[] k4 = multScalar(f4, delta_t);

        double[] res = new double[y_k.length];
        for (int i  = 0; i  < y_k.length; i ++)
            res[i] = y_k[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
        return res;
    }

    /**
     * addiert die zwei Vektoren a und b
     */
    private double[] addVectors(double[] a, double[] b) {
        double[] erg = new double[a.length];
        for (int i = 0; i < a.length; i++)
            erg[i] = a[i] + b[i];
        return erg;
    }

    /**
     * multipliziert den Skalar scalar auf den Vektor a
     */
    private double[] multScalar(double[] a, double scalar) {
        double[] erg = new double[a.length];
        for (int i = 0; i < a.length; i++)
            erg[i] = scalar * a[i];
        return erg;
    }

}
