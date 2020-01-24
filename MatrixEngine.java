public class MatrixEngine{
	
	public MatrixEngine() { }
	
	public static void UnitMatrix(double[][] A){ //Create unit matrix 
		int N = A.length;
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) { 
				if (i==j) A[i][j] = 1.0;
				else A[i][j] = 0.0;
	  		}
		}
	}
	
	public static void HilbertMatrix(double[][] A){ //Create Hilbert matrix 
		int N = A.length;
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) { 
				A[i][j]=1.0/(1.0+i+j);
	  		}
		} 
	}
	
	public static void RowColMatrix(double[][] A){ //Create RowCol matrix:
		// 1 2 3 4 5 6
		// 2 2 3 4 5 6
		// 3 3 3 4 5 6
		// 4 4 4 4 5 6
		// 5 5 5 5 5 6
		// 6 6 6 6 6 6
		int N = A.length;
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) { 
				if(i>j) A[i][j]=i+1;
				else A[i][j]=j+1;
			}
		}
	}
	
	public static void CoordinateMatrix(double[][] A){ //Create matrix each element of which
		// contains its own coordinates
		//for example, A[2][5] = 25
		int N = A.length;
		double Log10 = Math.log(N-1)/Math.log(10.0);
		int lg = (int)Log10; lg++;
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				A[i][j] = i*Math.pow(10,lg) +j;
			}
		}
	}
	
	public static void PrintMatrix(double [][] M){ //Print the square matrix content on console String S;
		int ROWs = M.length;
		int COLs = M[1].length;
		System.out.println("ROWs="+ROWs+" COLUMNs="+COLs); 
		for(int i=0; i<ROWs; i++){
			String S = "I="; S+= i;
			for(int j=0; j<COLs; j++){ 
				S+= " "+M[i][j];  	
			}	
			System.out.println(S); 
		} 
	}
	
	public static void TransposeMatrix(double[][] A){ //Interchange rows and columns of an initial cquare matrix
		int N = A.length;
		double W;
		for (int j=0; j<N-1; j++) {
			for (int i=j+1; i<N; i++) {
				W = A[i][j];
				A[i][j] = A[j][i];
				A[j][i] = W;
			}
		}
	}
	
	public static void TriangularMatrixUP(double[][] A){ //Create Triangular UP
		//matrix by zeroing
		//elements of 
		// an initial matrix:
		// ******
		// 0*****
		// 00****
		// 000***
		// 0000**
		// 00000*
		int N = A.length;
		for (int j=0; j<N-1; j++)
			for (int i=j+1; i<N; i++) 
				A[i][j]=0.0;
	}
	
	// SLAE creator
	public static void RightSum(double[][] A, double [][]B){ //Create Right Hand Parts as a SUM of Row
		// A - input matrix of size NxN;
		// B - output matrix of size NxN+1
		// B[0:N-1][0:N-1]=A[0:N-1][0:N-1]
		// B[I][N]=SUM(A[I][0...N-1])
		int N = A.length; 
		for (int i=0; i<N; i++)
			for (int j=0; j<N; j++)
				B[i][j]=A[i][j];
		
		for (int i=0; i<N; i++){
			B[i][N] = 0.0;
			for (int j=0; j<N; j++) 
				B[i][N]+=A[i][j];
		}
	}
	
	// Multiplier 01
	public static void MulMatrices(double[][] A, double[][] B, double[][] C){// Multiplication of two double  matrices 
		// C = A * B, where A and B are square 
		// initial matrices filled
		int Size = A.length;
		for (int i=0; i<Size; i++) {
			for (int j=0; j<Size; j++) { 
				C[i][j] = 0.0;
				for (int k=0; k<Size; k++)
					C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
	
	// Multiplier 02
	public static void MulMatricesMod(double[][] A, double[][] B,double[][] C){// Multiplication of two double matrices 
		// C = A * B, where A and B are square 
		// initial matrices filled
		// Modified Algorithm using middle-product method
		// By Row
		int Size = A.length;
		for (int i=0; i<Size; i++)
			for (int j=0; j<Size; j++)
				C[i][j] = 0.0;
		
		for (int i=0; i<Size; i++) {
			for (int j=0; j<Size; j++) {
				for (int k=0; k<Size; k++)
					C[j][k] += A[j][i]*B[i][k];
			}
		}
	}
	
	// Multiplier 03
	public static void MulMatricesRow(double[][] A, double[][] B, double[][] C){// Multiplication of two double matrices 
		// C = A * B, where A and B are square 
		// initial matrices filled. 
		int Size = A.length;
		TransposeMatrix(B);
		for (int i=0; i<Size; i++)
			for (int j=0; j<Size; j++){
				C[i][j] = 0.0;
				for (int k=0; k<Size; k++)
					C[i][j] += A[i][k]*B[j][k];
			}
		TransposeMatrix(B);
	} 
	
	// Solver of SLAE with Triangular Matrix
	public static void GetTriangularSolutionUP(double[][] A){// Solution of the matrix SLAE with the 
		// upper triangular matrix:
		// * * * * *
		// 0 * * * *
		// 0 0 * * *
		// 0 0 0 * *
		// 0 0 0 0 *
		int N = A.length; 
		int M = A[0].length;
		double S;
		for(int k=N; k<M; k++){ 
			A[N-1][k]=A[N-1][k]/A[N-1][N-1];
			for (int i=N-2; i>=0; i--) { 
				S=0;
				for (int j=i+1; j<N; j++){
					S+=A[i][j]*A[j][k];
				}
				A[i][k]=(A[i][k]-S)/A[i][i];
			}
		}
	}
	
	// Estimate error
	public static double GetErrorFromOnes(double[][] A){// Estimate error of the SLAE solving when
		// Solutions X[i] are ones;
		int N = A.length;
		int M = A[0].length;
		double ERROR=0.0;
		double ErrorMax = 0.0;
		for (int i=0; i<N; i++)
			for (int j=N; j<M; j++) { 
				ERROR = 1.0 - A[i][j];
				if (ERROR < 0)ERROR = - ERROR;
				if (ErrorMax < ERROR) ErrorMax = ERROR;
			}
		return ErrorMax;
	}
}