package alsCplex;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class ALS_NMF {

	private ArrayList<ArrayList<Integer> > neigh;
	
	private int n;
	private int m;
	private int k;
	private int NUMITER;
	private Random rand;
	private double lambda;
	private double[][] V;
	
	 /**
     * Creates an ALS_NMF object to stare the ALS_NMF algorithm 
     * @param n size of the input matrix (nxm)
     * @param m size of the input matrix (nxm)
     * @param k number of sources to separate
     * @param nIter number of iterations of the ALS algorithm
     * @param lambda regularization value for the neighborhood regularization
     * @param V input matrix
     * @param neigh list containing the "neighbor entries"
     * @param seed seed to initializate the first random matrix
     * 
     * @return an instance of ALS_NMF
     */
	
	public ALS_NMF(int n, int m, int k, int nIter, double lambda, double[][] V, ArrayList<ArrayList<Integer>> neigh, int seed) throws IOException{		
		this.n = n;
		this.m = m;
		this.V = V;
		this.neigh = neigh;
		this.k = k;
		NUMITER = nIter;
		//long randSeed = (long) (Math.random() * 100);
		rand = new Random(seed);
		this.lambda = lambda;
	}
	
	 /**
     * Executes the ALS_NMF algorithm 
     * 
     * @return an array list containing W and H matrices
     */
	
	public ArrayList<double[][]> executeNMF() throws IloException{
		ArrayList<double[][]> retVal = new ArrayList<double[][]>();
		double[][] wi = null;
		double[][] hi = null;
		
		hi = new double[k][];
		for(int i=0; i<k; i++)
			hi[i] = new double[m];
		
		for(int i=0; i<k; i++)
			for(int j=0; j<m; j++)
				hi[i][j] = rand.nextDouble();
		
		for(int iter=1; iter<=NUMITER; iter++){
			IloCplex cplex = new IloCplex();
			cplex.setOut(null);
			if(iter%2 == 1){
				IloNumVar[][] W = new IloNumVar[n][];
				for(int i=0; i<n; i++)
					W[i] = cplex.numVarArray(k, 0, Double.MAX_VALUE);
				wi = solveW(lambda, V, n, m, k, cplex, W, hi, neigh);
			}else{
				IloNumVar[][] H = new IloNumVar[k][];
				for(int i=0; i<k; i++)
					H[i] = cplex.numVarArray(m, 0, Double.MAX_VALUE);
				hi = solveH(lambda, V, n, m, k, cplex, H, wi, neigh);
			}
			cplex.end();
		}
		retVal.add(wi);
		retVal.add(hi);
		return retVal;
	}
	
	private static IloNumExpr[][] multiplyW(IloCplex cplex, IloNumVar[][] A, double[][] B) throws IloException {
		int mA = A.length;
		int nA = A[0].length;
		int mB = B.length;
		int nB = B[0].length;
		if (nA != mB) throw new RuntimeException("Illegal matrix dimensions.");
		IloNumExpr[][] C = new IloNumExpr[mA][];
		for(int i=0; i<mA; i++)
			C[i] = new IloNumExpr[nB];
		for (int i = 0; i < mA; i++){
			for (int j = 0; j < nB; j++){
				IloNumExpr expr = cplex.diff(A[0][0], A[0][0]);
				for (int w = 0; w < nA; w++)
					expr = cplex.sum(expr,cplex.prod(A[i][w],B[w][j]));
				C[i][j] = expr;
			}
		}
		return C;
	}
	
	private static IloNumExpr[][] multiplyH(IloCplex cplex, double[][] A, IloNumVar[][] B) throws IloException {
		int mA = A.length;
		int nA = A[0].length;
		int mB = B.length;
		int nB = B[0].length;
		if (nA != mB) throw new RuntimeException("Illegal matrix dimensions.");
		IloNumExpr[][] C = new IloNumExpr[mA][];
		for(int i=0; i<mA; i++)
			C[i] = new IloNumExpr[nB];
		for (int i = 0; i < mA; i++){
			for (int j = 0; j < nB; j++){
				IloNumExpr expr = cplex.diff(B[0][0], B[0][0]);
				for (int w = 0; w < nA; w++)
					expr = cplex.sum(expr,cplex.prod(A[i][w],B[w][j]));
				C[i][j] = expr;
			}
		}
		return C;
	}

	private static double[][] solveW(double lambda, double[][] V, int n, int m, int k, IloCplex cplex, IloNumVar[][] W, double[][] H, ArrayList<ArrayList<Integer> > neigh) throws IloException{
		IloNumExpr[][] WH = multiplyW(cplex, W, H);
		//Frobenius Norm
		IloNumExpr[] exp = new IloNumExpr[n*m];
		int count = 0;
		for(int i=0; i<n; i++){
			for(int j=0; j<m; j++){
				exp[count] = cplex.square(cplex.diff(V[i][j],WH[i][j]));
				count++;
			}
		}
		IloNumExpr objective = exp[0];
		for(int i=1; i<exp.length; i++)
			objective = cplex.sum(objective,exp[i]);

		//Spatial regularization
		if(lambda > 0){
			double regSum = 0;
			for(int i=0; i<neigh.size(); i++){
				for(int j=0; j<k; j++){
					int idx1 = neigh.get(i).get(0)-1;
					int idx2 = neigh.get(i).get(1)-1;
					regSum += (H[j][idx1]-H[j][idx2])*(H[j][idx1]-H[j][idx2]);
				}
			}
			regSum = (lambda/2)*regSum;
			objective = cplex.sum(objective,regSum);
		}
		//Objective function definition
		//System.out.println(objective);
		cplex.addMinimize(objective);

		//Options setting
		cplex.setParam(IloCplex.Param.Simplex.Display, 0);
		cplex.setParam(IloCplex.IntParam.NodeAlg, IloCplex.Algorithm.Dual);
		cplex.setParam(IloCplex.IntParam.SiftAlg, IloCplex.Algorithm.Dual);
		cplex.setDeleteMode(IloCplex.DeleteMode.LeaveBasis);
		cplex.setParam(IloCplex.Param.OptimalityTarget, 1);
		cplex.setParam(IloCplex.IntParam.PreDual, -1);
		cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Dual);

		//solve
		cplex.solve();
		//if(cplex.solve()) 
		//	System.out.println("obj = "+cplex.getObjValue());
		//else
		//	System.out.println("Model not solved");
		
		double[][] retVal = new double[n][];
		for(int i=0; i<n; i++)
			retVal[i] = new double[k];
		
		for(int i=0; i<n; i++)
			for(int j=0; j<k; j++)
				retVal[i][j] = cplex.getValue(W[i][j]);
		
		return retVal;
	}
	
	private static double[][] solveH(double lambda, double[][] V, int n, int m, int k, IloCplex cplex, IloNumVar[][] H, double[][] W, ArrayList<ArrayList<Integer> > neigh) throws IloException{
		IloNumExpr[][] WH = multiplyH(cplex, W, H);
		//Frobenius Norm
		IloNumExpr[] exp = new IloNumExpr[n*m];
		int count = 0;
		for(int i=0; i<n; i++){
			for(int j=0; j<m; j++){
				exp[count] = cplex.square(cplex.diff(V[i][j],WH[i][j]));
				count++;
			}
		}
		IloNumExpr objective = exp[0];
		for(int i=1; i<exp.length; i++)
			objective = cplex.sum(objective,exp[i]);

		
		//Spatial regularization
		IloNumExpr reg = cplex.diff(H[0][0], H[0][0]);
		for(int i=0; i<neigh.size(); i++){
			for(int j=0; j<k; j++){
				int idx1 = neigh.get(i).get(0)-1;
				int idx2 = neigh.get(i).get(1)-1;
				reg = cplex.sum(reg,cplex.square(cplex.diff(H[j][idx1], H[j][idx2])));
			}
		}
		reg = cplex.prod(lambda,reg);
		objective = cplex.sum(objective,reg);
		
		//Objective function definition
//		System.out.println(objective);
		cplex.addMinimize(objective);
		
		//Options setting
		cplex.setParam(IloCplex.Param.Simplex.Display, 0);
		cplex.setParam(IloCplex.IntParam.NodeAlg, IloCplex.Algorithm.Dual);
		cplex.setParam(IloCplex.IntParam.SiftAlg, IloCplex.Algorithm.Dual);
		cplex.setDeleteMode(IloCplex.DeleteMode.LeaveBasis);
		cplex.setParam(IloCplex.Param.OptimalityTarget, 1);
		cplex.setParam(IloCplex.IntParam.PreDual, -1);
		cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Dual);

		//solve
		cplex.solve();
		//if (cplex.solve()) 
		//	System.out.println("obj = "+cplex.getObjValue());
		//else
		//	System.out.println("Model not solved");
		
		double[][] retVal = new double[k][];
		for(int i=0; i<k; i++)
			retVal[i] = new double[m];
		
		for(int i=0; i<k; i++){
			for(int j=0; j<m; j++){
				retVal[i][j] = cplex.getValue(H[i][j]);
			}
		}
		return retVal;
	}
	
//	public static double[][] getBestSolution(ArrayList<ArrayList<double[][]>> inputResults){
//		int maxIdx = 0;
//		int maxSum = 0;
//		for(int i=0; i<inputResults.size(); i++){
//			double[][] matrix = inputResults.get(i).get(1);
//			int counter = countZero(matrix);
//			if(counter/((matrix.length*matrix[0].length)*2) < 0.5){
//				int tempSum = 0;
//				for(int j=0; j<matrix[0].length; j++)
//					tempSum += getMax(matrix, j)/getSum(matrix, j);
//				if(tempSum > maxSum){
//					maxSum = tempSum;
//					maxIdx = i;
//				}
//			}else{
//				System.out.println("Not accepted");
//			}
//		}	
//		return inputResults.get(maxIdx).get(1);
//	}
//	
//	private static double getMax(double[][] array, int columnIndex){
//		double max = 0;
//		for(int i=0; i<array.length; i++){
//			if(array[i][columnIndex] > max)
//				max = array[i][columnIndex];
//		}
//		return max;
//	}
//	
//	public static double getSum(double[][] array, int columnIndex){
//		double sum = 0;
//		for(int i=0; i<array.length; i++)
//			sum += array[i][columnIndex];
//		return sum;
//	}
//	
//	public static int countZero(double[][] matrix){
//		int count = 0;
//		for(int i=0; i<matrix.length; i++){
//			for(int j=0; j<matrix[0].length; j++){
//				if(matrix[i][j] == 0)
//					count++;
//			}
//		}
//		return count;
//	}
	
}
