package ode;

/**
 * Das Einschrittverfahren von Heun
 *
 * @author braeckle
 *
 */
public class Heun implements Einschrittverfahren {

    @Override
    /**
     * {@inheritDoc}
     * Nutzen Sie dabei geschickt den Expliziten Euler.
     */
    public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
        // TODO: done
        double[] result = new double[y_k.length];
        double[] f = ode.auswerten(t, y_k);

        ExpliziterEuler euler = new ExpliziterEuler();
        double[] euler_y_k = euler.nextStep(y_k,t,delta_t,ode);
        euler_y_k = ode.auswerten(t + delta_t, euler_y_k);

        for(int i = 0; i<result.length; i++)
            result[i] = y_k[i] + delta_t/2 * (euler_y_k[i] + f[i]);

        return result;
    }

}
