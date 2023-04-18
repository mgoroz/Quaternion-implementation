import java.lang.StringBuilder;
import java.lang.Override;
import java.util.Objects;


/** Quaternions. Basic operations. */
public class Quaternion {

   private double a;
   private double b;
   private double c;
   private double d;
   private static final double EPSILON = 1e-10;

   /** Constructor from four double values.
    * @param a real part
    * @param b imaginary part i
    * @param c imaginary part j
    * @param d imaginary part k
    */
   public Quaternion (double a, double b, double c, double d) {
      this.a = a;
      this.b = b;
      this.c = c;
      this.d = d;
   }




   /** Real part of the quaternion.
    * @return real part
    */
   public double getRpart() {
      return a;
   }

   /** Imaginary part i of the quaternion.
    * @return imaginary part i
    */
   public double getIpart() {
      return b;
   }

   /** Imaginary part j of the quaternion.
    * @return imaginary part j
    */
   public double getJpart() {
      return c;
   }

   /** Imaginary part k of the quaternion.
    * @return imaginary part k
    */
   public double getKpart() {
      return d;
   }

   /** Conversion of the quaternion to the string.
    * @return a string form of this quaternion:
    * "a+bi+cj+dk"
    * (without any brackets)
    */
   @Override
   public String toString() {
      StringBuilder sb = new StringBuilder();
      sb.append(a);
      if (b >= 0) {
         sb.append("+");
      }
      sb.append(b).append("i");
      if (c >= 0) {
         sb.append("+");
      }
      sb.append(c).append("j");
      if (d >= 0) {
         sb.append("+");
      }
      sb.append(d).append("k");
      return sb.toString();
   }



   /** Conversion from the string to the quaternion. 
    * Reverse to <code>toString</code> method.
    * @throws IllegalArgumentException if string s does not represent
    *     a quaternion (defined by the <code>toString</code> method)
    * @param s string of form produced by the <code>toString</code> method
    * @return a quaternion represented by string s
    */
   public static Quaternion valueOf(String s) {
      String[] parts = s.split("(?=[+-])"); // split string at + or - signs
      if (parts.length != 4) {
         throw new IllegalArgumentException("Invalid quaternion string: " + s);
      }
      double a, b, c, d;
      try {
         a = Double.parseDouble(parts[0]);
         b = Double.parseDouble(parts[1].substring(0, parts[1].length() - 1));
         c = Double.parseDouble(parts[2].substring(0, parts[2].length() - 1));
         d = Double.parseDouble(parts[3].substring(0, parts[3].length() - 1));
      } catch (NumberFormatException e) {
         throw new IllegalArgumentException("Invalid quaternion string: " + s, e);
      }
      return new Quaternion(a, b, c, d);
   }


   /** Clone of the quaternion.
    * @return independent clone of <code>this</code>
    */
   @Override
   public Object clone() throws CloneNotSupportedException {
      return new Quaternion(a, b, c, d);
   }



   /** Test whether the quaternion is zero. 
    * @return true, if the real part and all the imaginary parts are (close to) zero
    */
   public boolean isZero() {
      return a * a + b * b + c * c + d * d < EPSILON;
   }


   /** Conjugate of the quaternion. Expressed by the formula 
    *     conjugate(a+bi+cj+dk) = a-bi-cj-dk
    * @return conjugate of <code>this</code>
    */
   public Quaternion conjugate() {
      return new Quaternion(a, -b, -c, -d);
   }

   /** Opposite of the quaternion. Expressed by the formula 
    *    opposite(a+bi+cj+dk) = -a-bi-cj-dk
    * @return quaternion <code>-this</code>
    */
   public Quaternion opposite() {
      return new Quaternion(-a, -b, -c, -d);
   }

   /** Sum of quaternions. Expressed by the formula 
    *    (a1+b1i+c1j+d1k) + (a2+b2i+c2j+d2k) = (a1+a2) + (b1+b2)i + (c1+c2)j + (d1+d2)k
    * @param q addend
    * @return quaternion <code>this+q</code>
    */
   public Quaternion plus(Quaternion q) {
      return new Quaternion(a + q.a, b + q.b, c + q.c, d + q.d);
   }

   /** Product of quaternions. Expressed by the formula
    *  (a1+b1i+c1j+d1k) * (a2+b2i+c2j+d2k) = (a1a2-b1b2-c1c2-d1d2) + (a1b2+b1a2+c1d2-d1c2)i +
    *  (a1c2-b1d2+c1a2+d1b2)j + (a1d2+b1c2-c1b2+d1a2)k
    * @param q factor
    * @return quaternion <code>this*q</code>
    */
   public Quaternion times(Quaternion q) {
      double newA = a * q.a - b * q.b - c * q.c - d * q.d;
      double newB = a * q.b + b * q.a + c * q.d - d * q.c;
      double newC = a * q.c - b * q.d + c * q.a + d * q.b;
      double newD = a * q.d + b * q.c - c * q.b + d * q.a;
      return new Quaternion(newA, newB, newC, newD);
   }

   /** Multiplication by a coefficient.
    * @param r coefficient
    * @return quaternion <code>this*r</code>
    */
   public Quaternion times(double r) {
      return new Quaternion(a * r, b * r, c * r, d * r);
   }

   /** Inverse of the quaternion. Expressed by the formula
    *     1/(a+bi+cj+dk) = a/(a*a+b*b+c*c+d*d) + 
    *     ((-b)/(a*a+b*b+c*c+d*d))i + ((-c)/(a*a+b*b+c*c+d*d))j + ((-d)/(a*a+b*b+c*c+d*d))k
    * @return quaternion <code>1/this</code>
    */
   public Quaternion inverse() {
      double normSq = a * a + b * b + c * c + d * d;
      if (Math.abs(normSq) < EPSILON) {
         throw new ArithmeticException("Quaternion has zero norm and no inverse.");
      }
      double invNormSq = 1.0 / normSq;
      return new Quaternion(round(a * invNormSq, 4), round(-b * invNormSq, 4), round(-c * invNormSq, 4), round(-d * invNormSq, 4));
   }

   /** Difference of quaternions. Expressed as addition to the opposite.
    * @param q subtrahend
    * @return quaternion <code>this-q</code>
    */
   public Quaternion minus(Quaternion q) {
      return new Quaternion(a - q.a, b - q.b, c - q.c, d - q.d);
   }

   /** Right quotient of quaternions. Expressed as multiplication to the inverse.
    * @param q (right) divisor
    * @return quaternion <code>this*inverse(q)</code>
    */
   public Quaternion divideByRight (Quaternion q) {
      if (q.isZero()) {
         throw new RuntimeException("Cannot divide by zero!");
      }
      Quaternion conjQ = q.conjugate();
      Quaternion numerator = this.times(conjQ);
      double denominator = q.norm() * q.norm();
      Quaternion result = numerator.times(1.0 / denominator);
      return new Quaternion(
              result.getRpart(),
              result.getIpart(),
              result.getJpart(),
              result.getKpart()
      );

   }

   /** Left quotient of quaternions.
    * @param q (left) divisor
    * @return quaternion <code>inverse(q)*this</code>
    */
   public Quaternion divideByLeft (Quaternion q) {
      if (q.isZero()) {
         throw new RuntimeException("Cannot divide by zero!");
      }
      Quaternion conjQ = q.conjugate();
      Quaternion numerator = conjQ.times(this);
      double denominator = q.norm() * q.norm();
      Quaternion result = numerator.times(1.0 / denominator);
      return new Quaternion(
              result.getRpart(),
              result.getIpart(),
              result.getJpart(),
              result.getKpart()
      );

   }

   private double roundToDecimalPlaces(double d, int decimalPlaces) {
      double factor = Math.pow(2, decimalPlaces);
      return Math.round(d * factor) / factor;
   }

   private static double round(double value, int places) {
      if (places < 0) {
         throw new IllegalArgumentException("The number of decimal places must be non-negative.");
      }
      double scale = Math.pow(10, places);
      return Math.round(value * scale) / scale;
   }



   /** Equality test of quaternions. Difference of equal numbers
    *     is (close to) zero.
    * @param qo second quaternion
    * @return logical value of the expression <code>this.equals(qo)</code>
    */


   @Override
   public boolean equals(Object qo) {
      if (this == qo) {
         return true;
      }
      if (qo instanceof Quaternion) {
         Quaternion q = (Quaternion) qo;
         return Math.abs(q.a - a) < EPSILON && //need
                 Math.abs(q.b - b) < EPSILON &&
                 Math.abs(q.c - c) < EPSILON &&
                 Math.abs(q.d - d) < EPSILON;
      }
      return false;
   }




   /** Dot product of quaternions. (p*conjugate(q) + q*conjugate(p))/2
    * @param q factor
    * @return dot product of this and q
    */
   public Quaternion dotMult(Quaternion q) {
      Quaternion pConj = conjugate();
      Quaternion qConj = q.conjugate();
      Quaternion pqConj = pConj.times(q).plus(qConj.times(this)).times(0.5);
      return pqConj.conjugate();
   }

   /** Integer hashCode has to be the same for equal objects.
    * @return hashcode
    */
   @Override
   public int hashCode() {
      return Objects.hash(a, b, c, d);
   }



   /** Norm of the quaternion. Expressed by the formula 
    *     norm(a+bi+cj+dk) = Math.sqrt(a*a+b*b+c*c+d*d)
    * @return norm of <code>this</code> (norm is a real number)
    */
   public double norm() {
      return Math.sqrt(a * a + b * b + c * c + d * d);
   }


   public Quaternion pow (int n){
      if (n == 0) {
         return new Quaternion(1, 0, 0, 0);
      } else if (n == 1) {
         return new Quaternion(this.getRpart(), this.getIpart(), this.getJpart(), this.getKpart());
      } else if (n == -1) {
         return this.inverse();
      } else if  (n > 1){
         return this.times (this.pow(n - 1));
      } else {
         return this.pow(n * -1).inverse();
      }
   }

   /** Main method for testing purposes.
    * @param arg command line parameters
    */
   public static void main (String[] args) {

      // Creating some quaternions for testing
      Quaternion q1 = new Quaternion(3, 4, 2, 1);
      Quaternion q2 = new Quaternion(-2, 1, 5, -3);
      Quaternion q3 = new Quaternion(-3, -4, -2, -1);
      Quaternion q4 = new Quaternion(0, 0, 0, 0);
/*
      // Testing equals method
      System.out.println(q1.equals(q1)); // true
      System.out.println(q1.equals(q2)); // false
      System.out.println(q1.equals(q3)); // false

      // Testing clone method
      Quaternion q1Clone = null;
      try {
         q1Clone = (Quaternion) q1.clone();
      } catch (CloneNotSupportedException e) {
         e.printStackTrace();
      }
      System.out.println(q1Clone.equals(q1)); // true
      System.out.println(q1Clone == q1); // false

      // Testing toString method
      System.out.println(q1.toString()); // "3.0 + 4.0i + 2.0j + 1.0k"
      System.out.println(q2.toString()); // "-2.0 + 1.0i + 5.0j - 3.0k"
      System.out.println(q3.toString()); // "-3.0 - 4.0i - 2.0j - 1.0k"
      System.out.println(q4.toString()); // "0.0 + 0.0i + 0.0j + 0.0k"

      // Testing arithmetic methods
      Quaternion sum = q1.plus(q2);
      System.out.println(sum.toString());

      Quaternion diff = q1.minus(q2);
      System.out.println(diff.toString());

      Quaternion prod = q1.times(q2);
      System.out.println(prod.toString());

      Quaternion inverse = q1.inverse();
      System.out.println(inverse.toString());

      Quaternion divRight = q1.divideByRight(q2);
      System.out.println(divRight.toString());

      Quaternion divLeft = q2.divideByLeft(q1);
      System.out.println(divLeft.toString());
*/
      System.out.println(q1.pow(0).toString()); // "1.0 + 0.0i + 0.0j + 0.0k"
      System.out.println(q1.pow(1).toString()); // "3.0 + 4.0i + 2.0j + 1.0k" (equal to q1)
      System.out.println(q1.pow(-1).toString()); // inverse of q1
      System.out.println(q1.pow(5).toString()); // q1 * q1
      System.out.println(q1.pow(-2).toString()); // (q1 * q1).inverse()


   }
}

