import java.awt.PointerInfo;
import java.math.BigInteger;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Random;

// Edmond Mbadu 
//In this program we assumed that the eavesdropper Eve has captured 
//p, G, and NaG. It is left to compute Na using the baby step giant step
// algorithm.
public class Ellitpic_Curve_DiffieHellman {

	public static void main(String[] args) {

		// The generator point on the elliptic curve
		String G = "4 11";
		// The values that Alice sends to bob
		String NaG[] = new String[100];

		// The fixed prime p ( 7 digits)s
		BigInteger p = new BigInteger("1000253");
		// this array will hold a 100 values ( potential secret values from alice)
		BigInteger Na[] = new BigInteger[100];
		// the coefficient of x on the elliptic curve
		BigInteger a = new BigInteger("3");
		// the potential values of a
		String PotentialNa[] = new String[100];
//		System.out.println(logBabyStepGiantStep("3861 1242", G, p, a));

		for (int i = 0; i < NaG.length; i++) {
			Na[i] = getRandomBigInteger();
			NaG[i] = point_multiplication(G, Na[i], p, a);
			long start = System.nanoTime();

//			PotentialNa[i] = brute_froce(NaG[i], G, p, a);
			PotentialNa[i] = logBabyStepGiantStep(NaG[i], G, p, a);
			long endtime = System.nanoTime();
			System.out.println((endtime - start) / 1000000);
		}

	}

	/**
	 * See https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication for
	 * more info about elliptic curve point multiplication
	 * 
	 * @param G  the generator point on the curve
	 * @param Na is the scalar to multiply the point G by. namely Na.G
	 * @param p  is the fixed prime p
	 * @param a  is the coefficient a from the elliptic curve of the form Y^2= x^3 +
	 *           ax+b
	 * @return a string representation of Na.G ( where the x coordinates and y
	 *         coordinates are seperated by a space).
	 */
	public static String point_multiplication(String G, BigInteger Na, BigInteger p, BigInteger a) {
		String Gi = G;
		// if multiplied by zero which is the point at infinity
		// return infinity. Since zero is an absorbing element for multiplication,
		// the point at infinity (which is the zero for ECC) is also an absorbing
		// element
		if (Na.equals(BigInteger.ZERO)) {
			return "infinity";
		}
		// Na.G means G+G+..... G Na times
		for (BigInteger i = BigInteger.ONE; i.compareTo(Na) == -1; i = i.add(BigInteger.ONE)) {

			Gi = point_addition(G, Gi, p, a);

		}

		return Gi;
	}

	/**
	 * see https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication for
	 * more info about elliptic curve point addition
	 * 
	 * @param P a point on the elliptic curve
	 * @param Q is another point on the elliptic curve
	 * @param p is the fixed prime p
	 * @param a is the coefficient a from the elliptic curve of the form Y^2= x^3 +
	 *          ax+b
	 * @return P+Q ( as defined from point addition on elliptic curve).
	 */
	public static String point_addition(String P, String Q, BigInteger p, BigInteger a) {

		if (P.equalsIgnoreCase("infinity") || Q.equalsIgnoreCase("infinity")) {
			if (!P.equalsIgnoreCase("infinity") && Q.equalsIgnoreCase("infinity")) {
				return P;
			} else if (P.equalsIgnoreCase("infinity") && !Q.equalsIgnoreCase("infinity")) {
				return Q;
			}
			return "infinity";
		}
		String[] data = P.split(" ");
		String[] data1 = Q.split(" ");
		BigInteger Px = new BigInteger(data[0]);
		BigInteger Py = new BigInteger(data[1]);
		BigInteger Qx = new BigInteger(data1[0]);
		BigInteger Qy = new BigInteger(data1[1]);

		if (!P.equalsIgnoreCase(Q) && Px.equals(Qx)) {
			return "infinity";

		} else if (P.equalsIgnoreCase(Q) && Py.equals(BigInteger.ZERO) && Qy.equals(BigInteger.ZERO)) {
			return "infinity";
		} else if (!P.equalsIgnoreCase(Q)) {
			BigInteger m = ((Qy.subtract(Py)).multiply((Qx.subtract(Px)).modInverse(p))).mod(p);
			BigInteger xr = ((m.pow(2)).subtract(Px).subtract(Qx)).mod(p);
			BigInteger yr = ((m.multiply(Px.subtract(xr))).subtract(Py)).mod(p);

			return xr.toString() + " " + yr;
		}
		// last scenario
		BigInteger m = ((Px.pow(2).multiply(BigInteger.valueOf(3))).add(a))
				.multiply((BigInteger.valueOf(2).multiply(Py)).modInverse(p));
		BigInteger xr = ((m.pow(2)).subtract(Px).subtract(Qx)).mod(p);
		BigInteger yr = ((m.multiply(Px.subtract(xr))).subtract(Py)).mod(p);
		return xr.toString() + " " + yr;

	}

	public static String brute_froce(String NaG, String G, BigInteger p, BigInteger a) {

		String potential = "";
		for (BigInteger i = BigInteger.ONE; i.compareTo(p) == -1; i = i.add(BigInteger.ONE)) {
			potential = point_multiplication(G, i, p, a);
			if (potential.equalsIgnoreCase(NaG))
				return i.toString();
		}
		return "not found";

	}

	/**
	 * 
	 * @return a random Big Integer. The upper bound is the fixed prime p.
	 */
	public static BigInteger getRandomBigInteger() {
		Random rand = new Random();
		BigInteger upperLimit = new BigInteger("1000253");
		BigInteger result;
		do {
			result = new BigInteger(upperLimit.bitLength(), rand); // (2^4-1) = 15 is the maximum value
		} while (result.compareTo(upperLimit) >= 0); // exclusive of 13
		return result;
	}

	/**
	 * 
	 * @param point: the given point on the elliptic curve
	 * @return the inverse of the point; that is, its reflection according to the x
	 *         axis In other words, the inverse consists of switching the sign of
	 *         the y coordinate of the point
	 */

	public static String Inverse(String point) {
		if (point.equalsIgnoreCase("infinity")) {
			return "infinity";
		}

		String[] data = point.split(" ");
		BigInteger px = new BigInteger(data[0]);
		BigInteger py = new BigInteger(data[1]);
		// switching the sign of the y coordinate happens here
		py = py.multiply(BigInteger.valueOf(-1));
		return px.toString() + " " + py.toString();
	}

	/**
	 * 
	 * @param NaG the point alice sends to Bob
	 * @param G   the generator point
	 * @param p   the fixed prime p
	 * @param a   the coefficient of x in the equation of the elliptic curve
	 * @return the private value held only by Alice
	 */
	public static String logBabyStepGiantStep(String NaG, String G, BigInteger p, BigInteger a) {

		BigInteger m = sqrt(p).add(BigInteger.ONE);
		HashMap<String, BigInteger> h = new HashMap<>();
		String basePow = G;
		for (BigInteger j = BigInteger.valueOf(0); j.compareTo(m) < 0; j = j.add(BigInteger.ONE)) {

			if (j.equals(BigInteger.ZERO)) {
				h.put("infinity", j);
			} else {
				h.put(basePow, j);
				basePow = point_addition(basePow, G, p, a);
			}

		}
		String list2 = "";
		BigInteger target;
		for (BigInteger i = BigInteger.valueOf(0); i.compareTo(m) < 0; i = i.add(BigInteger.ONE)) {

			String in = Inverse(point_multiplication(G, i.multiply(m), p, a));
			list2 = point_addition(NaG, in, p, a);
			target = (BigInteger) h.get(list2);
			if (target != null)
				return (i.multiply(m).add(target)).mod(p).toString();

		}
		return BigInteger.valueOf(-1).toString();
	}

	public static BigInteger sqrt(BigInteger x) {
		BigInteger div = BigInteger.ZERO.setBit(x.bitLength() / 2);
		BigInteger div2 = div;
		// Loop until we hit the same value twice in a row, or wind
		// up alternating.
		for (;;) {
			BigInteger y = div.add(x.divide(div)).shiftRight(1);
			if (y.equals(div) || y.equals(div2))
				return y;
			div2 = div;
			div = y;
		}
	}
}