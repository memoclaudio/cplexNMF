import ilog.concert.IloException;
import utils.FileReaderWriter;

import java.io.IOException;
import java.util.ArrayList;

import alsCplex.ALS_NMF;

public class Main_NMF {

	public static void main(String[] args) throws IloException, IOException {
		
		int k=2;
		int ALSITERATIONS = 100;
		int REPETITIONS = 20;
		double lambda = 0.01;
//		double lambda = Double.parseDouble(args[0]);
//		int REPETITIONS = Integer.parseInt(args[1]);
		
		ArrayList<ArrayList<double[][]>> allSolutions = new ArrayList<ArrayList<double[][]>>();
		ArrayList<ArrayList<Integer>> neigh = FileReaderWriter.read("/Users/stamile/Desktop/NMF/trash/indexSet");
		ArrayList<ArrayList<Double>> v1 = FileReaderWriter.valuesReader("/Users/stamile/Desktop/NMF/trash/X1");
		
		int n = v1.size();
		int m = v1.get(0).size();
		
		//FIRST STEP
		double[][] V = new double[n][];
		for(int i=0; i<n; i++)
			V[i] = new double[m];
		for(int i=0; i<n; i++)
			for(int j=0; j<m; j++)
				V[i][j] = v1.get(i).get(j);

		for(int i=0; i<REPETITIONS; i++){
			System.out.println("Repetition: " + i);
			ALS_NMF nmf = new ALS_NMF(n, m, k, ALSITERATIONS, lambda, V, neigh, i);
			ArrayList<double[][]> resNMF = nmf.executeNMF();
			allSolutions.add(resNMF);
		}
		
	}

}
