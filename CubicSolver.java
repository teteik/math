import java.io.*;
import java.util.*;
import java.util.Locale;

public class CubicSolver {
    private static double eps;
    private static double delta;

    public static void main(String[] args) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader("input.txt"));
             PrintWriter out = new PrintWriter("output.txt")) {

            String line;
            int testNumber = 1;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                Scanner sc = new Scanner(line).useLocale(Locale.US);
                eps = sc.nextDouble();
                delta = sc.nextDouble();
                double a = sc.nextDouble();
                double b = sc.nextDouble();
                double c = sc.nextDouble();

                out.printf("Test %d:%n", testNumber);
                List<Double> roots = solveCubic(a, b, c);
                for (double root : roots) {
                    double fx = Math.abs(f(root, a, b, c));
                    out.printf("  %.10f  (|f(x)| = %.2e)%n", root, fx);
                }
                testNumber++;
            }
        }
    }

    static double f(double x, double a, double b, double c) {
        return x * x * x + a * x * x + b * x + c;
    }

    static List<Double> solveCubic(double a, double b, double c) {
        List<Double> roots = new ArrayList<>();
        double D = 4 * a * a - 12 * b;

        if (D <= 0) {
            double f0 = f(0, a, b, c);
            Interval interval;
            if (Math.abs(f0) < eps) {
                roots.add(0.0);
            } else if (f0 < 0) {
                interval = findFiniteIntervalRight(a, b, c, 0.0);
                roots.add(bisect(a, b, c, interval.left, interval.right));
            } else {
                interval = findFiniteIntervalLeft(a, b, c, 0.0);
                roots.add(bisect(a, b, c, interval.left, interval.right));
            }
        } else {
            double sqrtD = Math.sqrt(D);
            double alpha = (-2 * a - sqrtD) / 6.0;
            double beta = (-2 * a + sqrtD) / 6.0;

            double fAlpha = f(alpha, a, b, c);
            double fBeta = f(beta, a, b, c);

            if (fAlpha > eps && fBeta < -eps) {
                Interval leftInt = findFiniteIntervalLeft(a, b, c, alpha);
                roots.add(bisect(a, b, c, leftInt.left, leftInt.right));

                roots.add(bisect(a, b, c, alpha, beta));

                Interval rightInt = findFiniteIntervalRight(a, b, c, beta);
                roots.add(bisect(a, b, c, rightInt.left, rightInt.right));
            }
            else if (fAlpha > eps && fBeta > eps) {
                Interval leftInt = findFiniteIntervalLeft(a, b, c, alpha);
                roots.add(bisect(a, b, c, leftInt.left, leftInt.right));
            }
            else if (fAlpha < -eps && fBeta < -eps) {
                Interval rightInt = findFiniteIntervalRight(a, b, c, beta);
                roots.add(bisect(a, b, c, rightInt.left, rightInt.right));
            }
            else if (fAlpha > eps && Math.abs(fBeta) < eps) {
                roots.add(beta);
                Interval leftInt = findFiniteIntervalLeft(a, b, c, alpha);
                roots.add(bisect(a, b, c, leftInt.left, leftInt.right));
            }
            else if (Math.abs(fAlpha) < eps && fBeta < -eps) {
                roots.add(alpha);
                Interval rightInt = findFiniteIntervalRight(a, b, c, beta);
                roots.add(bisect(a, b, c, rightInt.left, rightInt.right));
            }
            else if (Math.abs(fAlpha) < eps && Math.abs(fBeta) < eps) {
                roots.add((alpha + beta) / 2.0);
            }
            else if (Math.abs(fAlpha) < eps && fBeta > eps) {
                roots.add(alpha);
            }
            else if (Math.abs(fBeta) < eps && fAlpha < -eps) {
                roots.add(beta);
            }
        }

        roots.sort(Double::compareTo);
        List<Double> unique = new ArrayList<>();
        for (double r : roots) {
            if (unique.isEmpty() || Math.abs(r - unique.get(unique.size() - 1)) > eps) {
                unique.add(r);
            }
        }
        return unique;
    }

    static class Interval {
        double left, right;
        Interval(double l, double r) {
            left = l;
            right = r;
        }
    }

    static Interval findFiniteIntervalRight(double a, double b, double c, double x0) {
        double x = x0;
        while (f(x, a, b, c) < 0) {
            x += delta;
        }
        return new Interval(x - delta, x);
    }

    static Interval findFiniteIntervalLeft(double a, double b, double c, double x0) {
        double x = x0;
        while (f(x, a, b, c) > 0) {
            x -= delta;
        }
        return new Interval(x, x + delta);
    }

    static double bisect(double a, double b, double c, double left, double right) {
        double fLeft = f(left, a, b, c);
        double fRight = f(right, a, b, c);

        if (Math.abs(fLeft) < eps) return left;
        if (Math.abs(fRight) < eps) return right;

        while (right - left > eps) {
            double mid = (left + right) / 2.0;
            double fMid = f(mid, a, b, c);

            if (Math.abs(fMid) < eps) {
                return mid;
            }

            if (fLeft * fMid < 0) {
                right = mid;
            } else {
                left = mid;
                fLeft = fMid;
            }
        }
        return (left + right) / 2.0;
    }
}
