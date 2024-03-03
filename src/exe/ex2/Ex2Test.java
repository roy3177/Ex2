package exe.ex2;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2023, Ariel University,
 *  * Ex2: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex2 - 
 * It contains few testing functions for the polynum functions as define in Ex2.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex2Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};

	@Test 
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex2.f(po1, 0);
		double fx1 = Ex2.f(po1, 1);
		double fx2 = Ex2.f(po1, 2);
		assertEquals(fx0, 2, Ex2.EPS);
		assertEquals(fx1, 4, Ex2.EPS);
		assertEquals(fx2, 6, Ex2.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex2.add(po1, po2);
		double f1x = Ex2.f(po1, x);
		double f2x = Ex2.f(po2, x);
		double f12x = Ex2.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex2.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex2.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex2.mul(po2, minus1);
		double[] p1 = Ex2.add(p12, pp2);
		assertTrue(Ex2.equals(p1, po1));

		//Another test:
		double[]p2= {1.0,3.0,5.0};
		double[]p3= {7.0,9.0,8.0};
		double[]result= {8.0,12.0,13.0};
		double[]ans=Ex2.add(p3, p2);
		assertArrayEquals(ans,result);
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex2.add(po1, po2);
		double[] p21 = Ex2.add(po2, po1);
		assertTrue(Ex2.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex2.add(po1, Ex2.ZERO);
		assertTrue(Ex2.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex2.mul(po1, Ex2.ZERO);
		assertTrue(Ex2.equals(p1, Ex2.ZERO));
		//Another test:
		double[]p2= {2.0,3.0,5.0};//5.0x^2+3.0x+2.0
		double[]p3= {5.0,3.0,1.0};//1.0x^2+3.0x+5.0
		double[]result= {10.0,21.0,36.0,18.0,5.0};//5.0x^4+18.0x^3+36.0x^2+21.0x+10.0
		double[]ans=Ex2.mul(p3, p2);
		assertArrayEquals(ans,result);
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex2.mul(po1, po2);
		double[] p21 = Ex2.mul(po2, po1);
		assertTrue(Ex2.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex2.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex2.f(po1, x);
			double f2x = Ex2.f(po2, x);
			double f12x = Ex2.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex2.EPS);
		}
	}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex2.derivative(p); // 2x + 6
		double[] dp2 = Ex2.derivative(dp1); // 2
		double[] dp3 = Ex2.derivative(dp2); // 0
		double[] dp4 = Ex2.derivative(dp3); // 0
		assertTrue(Ex2.equals(dp1, pt));
		assertTrue(Ex2.equals(Ex2.ZERO, dp3));
		assertTrue(Ex2.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex2.poly(p);
		double[] p1 = Ex2.getPolynomFromString(sp);
		double[] p2 = Ex2.getPolynomFromString(sp2);
		boolean isSame1 = Ex2.equals(p1, p);
		boolean isSame2 = Ex2.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex2.poly(p1));
		//Another test:

		String s = "5.0x^3 + 4.5x^2 - 3.5x + 6.0";
		double[] ans= {6.0, -3.5, 4.5, 5.0};
		double[]result=Ex2.getPolynomFromString(s);
		assertArrayEquals(ans,result);
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex2.ZERO, {1+Ex2.EPS/2}, {1,2}};
		double[][] xx = {{-2*Ex2.EPS}, {1+Ex2.EPS*1.2}, {1,2,Ex2.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex2.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex2.equals(d1[i], xx[i]));
		}
		//Test options
		//Option 1:
		double[] p1 = {1.0, 5.0, 12.0}; //12.0x^2+5.0x+1.0  
		//Equals.
		double[] p2 = {1.0, 5.0, 12.0};//12.0x^2+5.0x+1.0
		assertTrue(Ex2.equals(p1, p2));

		//Option 2:
		double[] p3 = {1.0, 6.0, 3.0};//3.0x^2+6.0x+1.0
		//Are not equals.
		double[] p4 = {1.0, 2.0, 4.0};//4.0x^2+2.0x+1.0
		assertFalse(Ex2.equals(p3, p4));

	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex2.sameValue(po1,po2, x1, x2, Ex2.EPS);
		double rs2 = Ex2.sameValue(po2,po1, x1, x2, Ex2.EPS);
		assertEquals(rs1,rs2,Ex2.EPS);
		//Another example of test:
		double[]p1= {3,6}; //6x+3
		double[]p2= {-2,5};//5x-2
		double hituch=Ex2.sameValue(p1, p2, -50, 50, Ex2.EPS); //Check the hituch of the polynoms.
		assertEquals(-5,hituch,Ex2.EPS);

	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=0, x2=4;
		double a1 = Ex2.area(po1, po2, x1, x2, 100);
		double a2 = Ex2.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2,Ex2.EPS);


		//Another test:
		double[] p1 = {1.0, 2.0, 3.0}; // 3.0x^2+2.0x+1.0
		double[] p2 = {1.0, 2.0, 3.0}; // 3.0x^2+2.0x+1.0
		double x3 = 0.0;
		double x4 = 1.0;
		int numberOfTrapezoid = 100;
		double result = Ex2.area(p1, p2, x1, x2, numberOfTrapezoid);
		assertEquals(0.0, result, Ex2.EPS);

	}
	@Test
	/**
	 * Tests the PolynomFromPoints function
	 */
	public void testPolynomFromPoints() {
		double[] xValues = {1.0, 2.0, 3.0};
		double[] yValues = {1.0, 4.0, 9.0};
		double[] ePolynom= Ex2.PolynomFromPoints(xValues, yValues);
		assertArrayEquals(new double[]{0.0,0.0,1.0},ePolynom);

	}
	@Test
	/**
	 * Tests the poly function
	 */
	public void testpoly() {

		double[]poly1= {1.0,2.0,3.0,4.0};
		String ans1="4.0x^3+3.0x^2+2.0x+1.0"; //Example number 1.
		String result1=Ex2.poly(poly1);
		assertEquals(ans1,result1);

		double[]poly2= {1.0,8.0,7.0,6.0};
		String ans2="6.0x^3+7.0x^2+8.0x+1.0"; //Example number 2.
		String result2=Ex2.poly(poly2);
		assertEquals(ans2,result2);


	}
	@Test
	/**
	 * Tests the derivative function
	 */
	public void testDerivative() {
		double[] p1 = new double[]{1.0, 2.0, 3.0}; // represents the polynomial x^2 + 2x + 3
		double[] ans1 = new double[]{2.0, 6.0}; // represents the derivative polynomial 2x + 6
		assertArrayEquals(ans1, Ex2.derivative(p1));

		double[] p2 = new double[]{4.0, 0.0, 0.0, 1.0}; // represents the polynomial x^3 + 4
		double[] ans2 = new double[]{0.0, 0.0, 3.0}; // represents the derivative polynomial 3x^2
		assertArrayEquals(ans2, Ex2.derivative(p2));
	}
	@Test
	/**
	 * Tests the length function
	 */
	public void testLength() {
		//Option 1:
		double[]poly=new double[5];
		double x1=0;
		double x2=20;
		int numberOfSegments=4;

		double ans=Ex2.length(poly, x1, x2, numberOfSegments);
		assertEquals(ans,20);
		//Option 2:
		double[]poly2=new double[4];
		double x3=2;
		double x4=24;
		int numberofSegments2=6;

		double ans2=Ex2.length(poly2, x3, x4, numberofSegments2);
		assertEquals(ans2,22);


	}
	@Test
	/**
	 * Tests the root_rec function:
	 */
	public void testRootRect() {
		//Test 1:
		double[] p = {0.0, 3.0, 1.0}; // p(x) = x^2 +3x
		double x1 = -4.0;
		double x2 = -2.0;
		double ans = -3.0; 
		double result = Ex2.root_rec(p, x1, x2, Ex2.EPS);
		assertEquals(ans, result);
		//Test 2:
		double[]p1= {-9.0,0.0,1.0}; //p(x)=x^2-9
		double x3=6.0;
		double x4=0.0;
		double ans1=3.0;
		double result1=Ex2.root_rec(p1, x3, x4, Ex2.EPS);
		assertEquals(ans1, result1);

	}
}



