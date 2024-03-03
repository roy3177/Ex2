package exe.ex2;

/**
 * Introduction to Computer Science 2023, Ariel University,
 * Ex2: arrays, static functions and JUnit
 *
 * This class represents a set of functions on a polynom - represented as array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynom: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code here")
 *
 * @author boaz.benmoshe
 */
public class Ex2 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynom is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};

	/**
	 * Computes the f(x) value of the polynom at x.
	 * @param poly
	 * @param x
	 * @return f(x) - the polynom value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans +=c*poly[i];
		}
		return ans;
	}
	/** Given a polynom (p), a range [x1,x2] and an epsilon eps. 
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x1) <= 0. 
	 * This function should be implemented recursively.
	 * @param p - the polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1); 
		double f2 = f(p,x2);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (f1*f2<=0 && Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
			//Calculating the vertex of a parabola given by three points:
			double mechane=((xx[0]-xx[1])*(xx[0]-xx[2])*(xx[1]-xx[2]));
			double a=((xx[0]*(yy[2]-yy[1]))+(xx[1]*(yy[0]-yy[2]))+(xx[2]*(yy[1]-yy[0])))/mechane;
			double b=((yy[1]-yy[0])/(xx[1]-xx[0]))-(a*(xx[0]+xx[1]));
			double c=yy[0]-(a*xx[0]*xx[0])-(b*xx[0]);

			double[]newArray=new double[3];

			newArray[0]=c;
			newArray[1]=b;
			newArray[2]=a;

			return newArray;
		}
		return ans;
	}
	/**
	 * Evaluates a polynomial with the given coefficients at the given value of x.
	 *This function works by Horner's method:
	 *The link for the method:https://en.wikipedia.org/wiki/Horner%27s_method.
	 * @param c the coefficients of the polynomial in descending order of degree
	 * @param x the value of x to evaluate the polynomial at
	 * @return the value of the polynomial at x
	 */
	public static double evaluate(double[]c,double x) {
		double ans = 0.0; //The  initial value of the polynomial at x.
		int n=c.length-1;
		for (int i = n; i >= 0; i--) {//Running in for loop in descending order of degree.
			ans = ans * x + c[i]; //Goes up every coefficient with the running of the loop.
		}
		return ans;
	}

	/** Two polynoms are equal if and only if the have the same values f(x) for 1+n values of x, 
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynom
	 * @param p2 second polynom
	 * @return true iff p1 represents the same polynom as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
		int n=Math.max(p1.length-1, p2.length-1); 
		for(int i=0;i<=n;i++) { //The loop is running on the maximum degree of the polynomials.
			double x = 1.0 + i; 			//y1=f(x1),y2=f(x2).
			double y1 = evaluate(p1, x);	//The function evaluates the polynomials at 1,2,...n+1
			double y2 = evaluate(p2, x);	//Connection with the function evaluate.
			if (Math.abs(y1 - y2) > EPS) { 	
				return false;
			}
		}
		return ans;
	}

	/** 
	 * Computes a String representing the polynom.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynom represented as an array of doubles
	 * @return String representing the polynom: 
	 */
	public static String poly(double[] poly) {
		String ans = ""; 
		if(poly.length==0 || poly.length==1) {
			if(poly.length==0) {
				ans="0";	
			}
			else {
				ans=Double.toString(poly[0]);
			}
		}
		else {
			int l=poly.length-1; 
			for(int p=l;p>=0;p--) { //Running with the loop for in descending order
				if(Math.abs(poly[p])<EPS) {
					continue;				//p meaning the power of the numbers in the array.
				}
				if(ans.length()>0 && poly[p]>0) {
					ans=ans+ "+";				//connect between the chars in the String.
				}
				if(Math.abs(poly[p]-1)>EPS || p==0) {
					ans=ans+poly[p];              // Checking the coefficients.
				}
				if(p>0) {
					ans=ans+ "x";
				}
				if(p>1) {
					ans=ans+"^"+p; //There is a number in the power that is not 
				}				   //0,or 1,which should be expressed in a String.
			}


		}
		return ans;
	}
	/**
	 * Given two polynoms (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {

		// add you code here
		double Xm = x1;		 //The initial value of Xm(x1<=Xm<=x2).
		double y1 = f(p1, x1) - f(p2, x1);		 //The value of x1 between the polynoms.
		double y2 = f(p1, x2) - f(p2, x2);		//The value of x2 between the polynoms.
		if (y1 * y2 > 0) {            		//assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2))>0.
			return Double.NaN;     			   //Return invalid value 
		}
		while (Math.abs(x2 - x1) > eps) {
			Xm = (x1 + x2) / 2;		//The middle between of the two points of x:(x1,x2).
			double Ym = f(p1, Xm) - f(p2, Xm);	//The value of f(Xm) between p1 and p2.
			if (Math.abs(Ym) < eps) {	//|p1(x) -p2(x)| < epsilon.
				return Xm;				//Return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < epsilon.
			}
			if (Ym * y1 < 0) {
				x2 = Xm;		//The p2 got the p(Xm,Ym).
				y2 = Ym;
			} else {
				x1 = Xm;		//The p1 got the p(Xm,Ym).
				y1 = Ym;
			}
		}
		return Xm;
	}
	/**
	 * Given a polynom (p), a range [x1,x2] and an integer with the number (n) of sample points. 
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = 0.0; //Initial value of the length.

		double segmentWidth=(x2-x1)/(numberOfSegments); //Calculates the size of each segment.
		double xNext=segmentWidth+x1; //Calculate the next x in the length.

		for(double i=x1;i<x2;i=i+segmentWidth) {	//Running in the range (x1,x2).

			// Calculate the length of the segment using distance formula:
			double distance=Math.sqrt(Math.pow((xNext-i),2)+Math.pow((Ex2.f(p, xNext)-Ex2.f(p, i)), 2));
			ans=ans+distance;  // Add the segment length to the total length,from any distance.
			xNext=xNext+segmentWidth;	//Any iteration,we update the value of xNext.
		}

		return ans; //Return the final value of the length.
	}

	/**
	 * Given two polynoms (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom). 
	 * This function computes an approximation of the area between the polynoms within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynoms within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		if(numberOfTrapezoid<=0) {
			return 0;		//Cannot calculate the reimann's integral.
		}
		double ans = 0; //The initial value of the area.
		double width=Math.abs((x2-x1)/numberOfTrapezoid); //Calculate the width of the trapezoid.

		for(double i=x1;i<x2;i=i+width) {
			double inter=sameValue(p1,p2,i,i+width,EPS); //Checks the exciting intersection.
			if(inter>i && inter<i+width) {
				double y1=Math.abs(inter-i);		//Calculate the intersection points.
				double y2=Math.abs(i+width-inter);
				double l=f(p1,i)-f(p2,i);
				double abs=Math.abs(l);//Calculate the length of the boxes via absolute number for the functions.
				double l2=Math.abs(f(p1,i+width)-f(p2,i+width)); //Moves over the integral.
				ans=ans+(abs*y1*0.5)+(l2*y2*0.5);
			}
			else {
				double result=((Math.abs(f(p1,i)- f(p2,i))+Math.abs(f(p1,i+width)- f(p2,i+width)))*width)/2;
				ans=ans+result;
			}
		}
		return ans; //Return the final value of the area.
	}
	/**
	 * This function computes the array representation of a polynom from a String
	 * representation. Note:given a polynom represented as a double array,  
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynom.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		p = p.toLowerCase();
		String []s1= p.split("((?>=-)|(?=-)|(?>=\\+)|(?=\\+))");//Using the split method on the array.                               
		if(p.length()==0) {          //If the string empty.
			return new double[]{0.0}; 
		}
		double  number;
		String []r;
		int maxPower=0; //The initial value of the maximum power.
		for (int i=0; i<s1.length; i++) {                                                         
			s1[i]=s1[i].replaceAll(" ","");   //Delete the places the are empty.                                               
		}
		for (int i=0; i<s1.length; i++) {                           
			int l=s1[i].length(); //The new length after the first split.                                                       
			if(s1[i].contains("^")) { //Checking if there is a power in the string.                                                           
				char saverPower= s1[i].charAt(l-1);                                            
				int newPower= Integer.parseInt(String.valueOf(saverPower));   
				if(newPower>maxPower){                                                                
					maxPower=newPower;                                                               
				}
			}
		}	
		double[]a= new double [maxPower+1];    //Define array at the size of maxPower.                                                           
		for (int i=0; i<s1.length; i++) {                            
			if(s1[i].contains("^")) {                                                             
				r=s1[i].split("x"); //Removing x from the string.                                                   
				r[1]=r[1].substring(1);  // remove ^ from string                                                 
				double powerD= Double.parseDouble(r[1]);                                               
				int powerI= (int) powerD;   //Casting the power.                                                             
				number= Double.parseDouble(r[0]);                                              
				if (a[powerI]==0) {                                                 	           
					a[powerI]=number;  // Putting  the number in the empty cell.
				}
				else {
					a[powerI]=a[powerI]+number;  // Adding a number and giving new value in the cell.                                                    
				}
			}
			// If there is not a power.
			else {                                                                                 
				if (s1[i].contains("x")) {                                                        
					r=s1[i].split("x");                                                      
					number= Double.parseDouble(r[0]);                                            
					if (a[1]==0) {                                                               
						a[1]=number;                                                              
					}
					else {
						a[1]=a[1]+number;                                               
					}
				}
				//If there is only a number.
				else {                                                                             
					number=Double.parseDouble(s1[i]);	                                            
					if (a[0]==0) {                                                              
						a[0]=number;                                                                
					}
					else {            
						a[0]=a[0]+number;                                                          
					}                                                                
				}
			}
		}
		if (s1.length==0) {  //If the string is empty.
			a=null;                                                                              
		}		
		return a; 
	}

	/**													
	 * This function computes the polynom which is the sum of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
		// add you code here	
		int p=Math.max(p1.length, p2.length); //Findings the maximal length of p1 and p2(one of them).
		double[]sumPolynoms=new double[p];		

		if(p>0) {
			for(int i=0;i<p;i++) {
				if (i < p1.length) {
					sumPolynoms[i]=sumPolynoms[i]+p1[i];
				}
				if (i < p2.length) {
					sumPolynoms[i]=sumPolynoms[i]+ p2[i];
				}
			}
			return sumPolynoms; //Return the sum of the two polynoms (p1,p2).

		}
		else {
			return ans; //If two of the polynoms equals to zero.
		}
	}
	/**
	 * This function computes the polynom which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;
		int lp1 = p1.length; 	//The length of the polynoms.
		int lp2 = p2.length;

		if(lp2==0 || lp1==0) {
			return ans;			//Return empty array that equal to the ans array.
		}
		ans=new double[lp1+lp2-1];	//Give the ans array a new value with the lengths of the polynoms.
		for (int i = 0; i < lp1; i++) {		//Running on p1 and p2.
			for (int j = 0; j < lp2; j++) {			
				ans[i + j]=ans[i + j]+ p1[i] * p2[j];	// Multiply each place of p1 with each place of p2.
			}
		}
		return ans;
	}

	/**
	 * This function computes the derivative polynom:.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;

		int lp=po.length-1;		//  The degree of the input polynom minus one. 
		double[]d=new double[lp]; // The derivative polynom,comes in new array.
		if(lp!=0){
			for (int i = 1; i <= lp; i++) {	//Going over the coefficients in the array po.
				d[i - 1] = i * po[i];		//Compute the derivative.
			}
			return d;
		}
		return ans;
	}
}
