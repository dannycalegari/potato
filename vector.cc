/*	vector.cc basic linear algebra functions	*/


dvec INV(dmat M, dvec L){	// finds vector U[0] st. M*U=L. If M is invertible, this is M^{-1}L.

	/* This function finds a vector u such that M*u = L. It does this by starting with the
	identity matrix U and a copy N of M, and then adjusts N by column operations until it
	is equal to the  vector L. We perform the same sequence of column operations on U.
	Then the desired vector is u = U[0].	*/

	// WARNING: we assume M has full rank! //
	
	// M[i][j] means row i, column j
	
//	cout << "Inverting matrix. Step 0: initialization.\n";
	int rows, cols;
	int i,j,k,l;
	rows=(int) M.size();
	cols=(int) M[0].size();		// M has size rows x cols; L should be 1 x cols
	dmat N,U;
	dvec K;
	dbl w;
	dbl t;
	N=M;
	U.resize(0);
	for(i=0;i<cols;i++){		// U initialized to cols x cols identity matrix
		K.resize(0);
		for(j=0;j<cols;j++){
			if(i==j){
				K.push_back(1.0);
			} else {
				K.push_back(0.0);
			};
		};
		U.push_back(K);
	};
	K.resize(0);
	
	/* Step 1. Make N lower diagonal. For each i, first swap columns >= i 
		so that column i has the biggest entry in row i, then subtract multiples
		of this column from columns > i so that their entry in row i is zero.
		
	see J. von Neumann, H.H. Goldstine, Numerical Inverting of Matrices of High Order, 
        Bull. Amer. Math. Soc., Vol. 53, No. 11, pp. 1021-1099, 1947.	*/
        
//	cout << "Inverting matrix. Step 1: Triangularization.\n";

	for(i=0;i<rows;i++){
		k=i;
		t=abs(N[i][k]);
		for(j=i;j<cols;j++){		// find j>=i maximizing N[i][j]
			if(abs(N[i][j])>t){
				k=j;
				t=abs(N[i][j]);
			};
		};
		for(j=0;j<rows;j++){	// swap columns i and k of N
			w=N[j][i];
			N[j][i]=N[j][k];
			N[j][k]=w;
		};
		for(j=0;j<cols;j++){	// keep U in sync
			w=U[i][j];
			U[i][j]=U[k][j];
			U[k][j]=w;
		};		
		for(j=i+1;j<cols;j++){
			w=N[i][j]/N[i][i];
			for(l=0;l<rows;l++){	// add multiple of column i to column j>i to make N[i][j]=0
				N[l][j]=N[l][j]-w*N[l][i];
			};
			for(l=0;l<cols;l++){	// keep U in sync
				U[j][l]=U[j][l]-w*U[i][l];
			};
		};
	};
//	cout << "Inverting matrix. Step 2: Subtraction.\n";

	/* Step 2. Adjust the first column of N to be equal to L.	*/
	for(i=0;i<rows;i++){
		w=(L[i]-N[i][0])/N[i][i];
		for(j=0;j<rows;j++){	// add w times column i to column 0
			N[j][0]=N[j][0]+w*N[j][i];
		};
		for(j=0;j<cols;j++){	// keep U in sync
			U[0][j]=U[0][j]+w*U[i][j];
		};
	};
	
	return(U[0]);
};