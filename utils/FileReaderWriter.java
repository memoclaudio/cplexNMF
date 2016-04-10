package utils;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class FileReaderWriter {
	
	public static ArrayList<ArrayList<Integer> > read(String filename) throws IOException{
		ArrayList<ArrayList<Integer> > retVal = new ArrayList<ArrayList<Integer> >();
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line = null;
		while((line = br.readLine()) != null){
			ArrayList<Integer> tmp = new ArrayList<Integer>();
			String[] field = line.split(" ");
			for(int i=0; i<field.length; i++){
				int current = Integer.parseInt(field[i]);
				tmp.add(current);
			}
			retVal.add(tmp);
		}
		br.close();
		return retVal;
	}
	
	public static ArrayList<ArrayList<Double> > valuesReader(String filename) throws IOException{
		ArrayList<ArrayList<Double> > retVal = new ArrayList<ArrayList<Double> >();
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line = null;
		while((line = br.readLine()) != null){
			ArrayList<Double> tmp = new ArrayList<Double>();
			String[] field = line.split(" ");
			for(int i=0; i<field.length; i++){
				double current = Double.parseDouble(field[i]);
				tmp.add(current);
			}
			retVal.add(tmp);
		}
		br.close();
		return retVal;
	}

	public static void writeFile(String filename, IloCplex cplex, IloNumVar[][] X) throws FileNotFoundException, UnknownObjectException, IloException{
		int n = X.length;
		int m = X[0].length;
		PrintWriter writer = new PrintWriter(filename);
		
		for(int i=0; i<n; i++){
			for(int j=0; j<m; j++){
				if(j != m-1)
					writer.print(cplex.getValue(X[i][j]) + " ");
				else
					writer.print(cplex.getValue(X[i][j]));
			}
			writer.println();
		}
		writer.close();
	}
	
	public static void writeFile(String filename, double[][] X) throws FileNotFoundException{
		int n = X.length;
		int m = X[0].length;
		PrintWriter writer = new PrintWriter(filename);
		
		for(int i=0; i<n; i++){
			for(int j=0; j<m; j++){
				if(j != m-1)
					writer.print(X[i][j] + " ");
				else
					writer.print(X[i][j]);
			}
			writer.println();
		}
		writer.close();
	}
	
	public static void writeCoordinate(String filename, ArrayList<ArrayList<Integer>> X) throws FileNotFoundException{
		PrintWriter writer = new PrintWriter(filename);
		
		for(ArrayList<Integer> currCorrd : X)
			writer.println(currCorrd.get(0) + " " + currCorrd.get(1) + " " + currCorrd.get(2));
		
		writer.close();
	}
	
	public static double[][] readInitMatrix(String filename) throws NumberFormatException, IOException{
		ArrayList<ArrayList<Double> > retVal = new ArrayList<ArrayList<Double> >();
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line = null;
		while((line = br.readLine()) != null){
			ArrayList<Double> tmp = new ArrayList<Double>();
			String[] field = line.split(" ");
			for(int i=0; i<field.length; i++){
				double current = Double.parseDouble(field[i]);
				tmp.add(current);
			}
			retVal.add(tmp);
		}
		br.close();
		double[][] finalMatrix = new double[retVal.size()][];
		for(int i=0; i<retVal.size(); i++)
			finalMatrix[i] = new double[retVal.get(i).size()];
		for(int i=0; i<retVal.size(); i++)
			for(int j=0; j<retVal.get(i).size(); j++)
				finalMatrix[i][j] = retVal.get(i).get(j);
		return finalMatrix;
	}
}
