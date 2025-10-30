import java.util.function.DoubleUnaryOperator;
import java.util.*;

public class NumericalIntegration {
    public static double leftRectangle(DoubleUnaryOperator f, double a, double b, int n) {
        double h = (b - a) / n;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += f.applyAsDouble(a + i * h);
        }
        return sum * h;
    }

    public static double rightRectangle(DoubleUnaryOperator f, double a, double b, int n) {
        double h = (b - a) / n;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += f.applyAsDouble(a + (i + 1) * h);
        }
        return sum * h;
    }

    public static double midpoint(DoubleUnaryOperator f, double a, double b, int n) {
        double h = (b - a) / n;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += f.applyAsDouble(a + (i + 0.5) * h);
        }
        return sum * h;
    }

    public static double trapezoidal(DoubleUnaryOperator f, double a, double b, int n) {
        double h = (b - a) / n;
        double sum = f.applyAsDouble(a) + f.applyAsDouble(b);
        for (int i = 1; i < n; i++) {
            sum += 2.0 * f.applyAsDouble(a + i * h);
        }
        return sum * h / 2.0;
    }

    public static double simpson(DoubleUnaryOperator f, double a, double b, int n) {
        if (n % 2 == 1) n++; // Симпсон требует чётное n
        double h = (b - a) / n;
        double sum = f.applyAsDouble(a) + f.applyAsDouble(b);
        for (int i = 1; i < n; i++) {
            double coeff = (i % 2 == 0) ? 2.0 : 4.0;
            sum += coeff * f.applyAsDouble(a + i * h);
        }
        return sum * h / 3.0;
    }


    static class IntegrationMethod {
        String name;
        java.util.function.Function<Params, Double> func;

        IntegrationMethod(String name, java.util.function.Function<Params, Double> func) {
            this.name = name;
            this.func = func;
        }
    }

    static class Params {
        DoubleUnaryOperator f;
        double a, b;
        int n;
        Params(DoubleUnaryOperator f, double a, double b, int n) {
            this.f = f; this.a = a; this.b = b; this.n = n;
        }
    }


    public static void printTable(DoubleUnaryOperator f, double a, double b, String expr, int[] Ns) {
        System.out.println("\n=== " + expr + " on [" + a + ", " + b + "] ===");

        System.out.printf("%-12s", "Method");
        for (int n : Ns) {
            System.out.printf("%20s", "n=" + n);
        }
        System.out.println();

        System.out.println("-".repeat(12 + 20 * Ns.length));

        IntegrationMethod[] methods = {
                new IntegrationMethod("Left Rect", p -> leftRectangle(p.f, p.a, p.b, p.n)),
                new IntegrationMethod("Right Rect", p -> rightRectangle(p.f, p.a, p.b, p.n)),
                new IntegrationMethod("Midpoint", p -> midpoint(p.f, p.a, p.b, p.n)),
                new IntegrationMethod("Trapezoidal", p -> trapezoidal(p.f, p.a, p.b, p.n)),
                new IntegrationMethod("Simpson", p -> simpson(p.f, p.a, p.b, p.n))
        };

        for (IntegrationMethod m : methods) {
            System.out.printf("%-12s", m.name);
            for (int n : Ns) {
                try {
                    double val = m.func.apply(new Params(f, a, b, n));
                    System.out.printf("%20.10f", val);
                } catch (Exception e) {
                    System.out.printf("%20s", "error");
                }
            }
            System.out.println();
        }
    }


    public static double rungeOrder(java.util.function.Function<Params, Double> method,
                                    DoubleUnaryOperator f, double a, double b, int n0) {
        double I1 = method.apply(new Params(f, a, b, n0));
        double I2 = method.apply(new Params(f, a, b, 2 * n0));
        double I3 = method.apply(new Params(f, a, b, 4 * n0));
        double denom = I2 - I3;
        if (Math.abs(denom) < 1e-15) return Double.POSITIVE_INFINITY;
        return Math.log(Math.abs(I1 - I2) / Math.abs(denom)) / Math.log(2.0);
    }

    public static void testRunge() {
        DoubleUnaryOperator f1 = x -> Math.pow(x, 11) - 5 * Math.pow(x, 8) + 0.654 * Math.pow(x, 3) - x;
        DoubleUnaryOperator f2 = x -> Math.sin(x * x - 0.3 * x + 4 * Math.PI) *
                Math.cos(Math.PI / 7.0 * x + 1);

        double[][] intervals = {{-3, 1}, {-Math.PI, 15 * Math.PI / 2}};
        String[] names = {
                "x^11 - 5*x^8 + 0.654*x^3 - x",
                "sin(x^2 - 0.3*x + 4*Pi) * cos(Pi/7*x + 1)"
        };

        IntegrationMethod[] methods = {
                new IntegrationMethod("Left Rect", p -> leftRectangle(p.f, p.a, p.b, p.n)),
                new IntegrationMethod("Right Rect", p -> rightRectangle(p.f, p.a, p.b, p.n)),
                new IntegrationMethod("Midpoint", p -> midpoint(p.f, p.a, p.b, p.n)),
                new IntegrationMethod("Trapezoidal", p -> trapezoidal(p.f, p.a, p.b, p.n)),
                new IntegrationMethod("Simpson", p -> simpson(p.f, p.a, p.b, p.n))
        };

        int n0 = 10;

        for (int idx = 0; idx < 2; idx++) {
            DoubleUnaryOperator f = (idx == 0) ? f1 : f2;
            double a = intervals[idx][0], b = intervals[idx][1];
            System.out.println("\n=== Runge order estimation for " + names[idx] +
                    " on [" + a + ", " + b + "] ===");
            System.out.printf("%-12s%12s%n", "Method", "Order (p)");
            System.out.println("-".repeat(25));
            for (IntegrationMethod m : methods) {
                try {
                    double p = rungeOrder(m.func, f, a, b, n0);
                    if (Double.isInfinite(p)) {
                        System.out.printf("%-12s%12s%n", m.name, "inf");
                    } else {
                        System.out.printf("%-12s%12.4f%n", m.name, p);
                    }
                } catch (Exception e) {
                    System.out.printf("%-12s%12s%n", m.name, "error");
                }
            }
        }
    }


    public static void main(String[] args) {
        List<DoubleUnaryOperator> funcs = Arrays.asList(
                x -> 64 * x + 19.82,
                x -> Math.pow(x, 15) - 4 * Math.pow(x, 8) + Math.pow(x, 3) - 147,
                x -> Math.sqrt(x + 32),
                x -> 1.0 / Math.pow(x + 3, 5),
                x -> Math.log(x * x + 9 * x - 3),
                x -> Math.sin((6.0 / 7.0) * Math.PI * x + x - 3 * Math.PI),
                x -> 1.0 / (1.0 + Math.exp(x))
        );

        double[][] bounds = {
                {-2, 4},
                {-6, -2},
                {1, 2},
                {0, 3},
                {10, 12},
                {-Math.PI / 4, 9 * Math.PI / 5},
                {0, 1}
        };

        String[] exprs = {
                "64*x + 19.82",
                "x^15 - 4*x^8 + x^3 - 147",
                "sqrt(x + 32)",
                "1/(x + 3)^5",
                "ln(x^2 + 9*x - 3)",
                "sin(6/7*Pi*x + x - 3*Pi)",
                "1/(1 + exp(x))"
        };

        int[] Ns = {10, 20, 50, 100, 1000};

        for (int i = 0; i < funcs.size(); i++) {
            printTable(funcs.get(i), bounds[i][0], bounds[i][1], exprs[i], Ns);
        }

        testRunge();
    }
}
