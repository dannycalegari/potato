/* dual_format.cc reads data in dual (i.e. triangle) format and converts */

/* Here is a description of "dual format". A file is in the form:

VERTICES n
infty
FACES m
i j k l_jk l_ik l_ij
... (repeated m times)

meaning: there are n vertices; vertex infty (in {1...n}) is to be put at infinity; 
there are m faces, all triangular; each line describes a triangle. the line above 
describes the triangle with vertices i j k on its CCW boundary, and lengths of 
edges jk, ik, ij (not log-lengths!)	*/


void packing::dual_format_input_data(ifstream &input_file){
	string S;
	int vertices,faces,infinity;
	ivl FACE;	// vector of faces
	dvl LENGTH;	//
	ivec F;
	dvec FL;
	int i,j,k;
	int v;
	dbl l;
	clock_t start;
	start=clock();
	
		// Step 1: read data from file
		
	input_file >> S;
	assert(S=="VERTICES");
	input_file >> vertices;
	input_file >> infinity;
	input_file >> S;
	assert(S=="FACES");
	input_file >> faces;
	FACE.resize(0);
	LENGTH.resize(0);
	for(i=0;i<faces;i++){
		F.resize(0);
		FL.resize(0);
		for(j=0;j<3;j++){
			input_file >> v;
			F.push_back(v);
		};
		FACE.push_back(F);
		for(j=0;j<3;j++){
			input_file >> l;
			FL.push_back(l);
		};
		LENGTH.push_back(FL);
	};

		// Step 2: convert to standard format
		
	ivec I;	// scratch vector for ADJ 
	dvec T;	// scratch vector for L
	int m,next,last;
		
	ADJ.resize(0);	// initialize data
	L.resize(0);
	SWAP_LIST.resize(0);

	cout << "computing vertex adjacency lists \n";
	for(i=0;i<vertices;i++){
		// generate adjacency list for vertex i
		// NOTE!!!! in data_list vertices start at 1; in standard format they start at 0.
		// So vertex i in standard format corresponds to vertex i+1 is data_list format.
		I.resize(0);	// initialize scratch vectors
		T.resize(0);
		
		// find first edge
		for(j=0;j<faces;j++){
			for(k=0;k<3;k++){
				if(FACE[j][k]==i+1){	// did we find vertex i?
					m=FACE[j][(k+1)%3]-1;	// what do we point to?
					last=m;	
					next=FACE[j][(k+2)%3]-1;	// what is the next one to look for?
					I.push_back(m);		// this is the first edge
					T.push_back(2.0*log(LENGTH[j][(k+2)%3]));	// L variable of the edge
					k=3;		//	exit from loop
					j=faces;	//	exit from loop
				};
			};
		};
		
		// find subsequent edges
		while(next!=last){
			for(j=0;j<faces;j++){
				for(k=0;k<3;k++){
					if(FACE[j][k]==i+1 && FACE[j][(k+1)%3]==next+1){	// is this the next edge?
						m=FACE[j][(k+1)%3]-1;
						next=FACE[j][(k+2)%3]-1;	// adjust value of next
						I.push_back(m);		// this is the first edge
						T.push_back(2.0*log(LENGTH[j][(k+2)%3]));	// L variable of the edge
						k=3;		//	exit from loop
						j=faces;	//	exit from loop								
					};
				};
			};
		};		
		ADJ.push_back(I);
		L.push_back(T);
	};
	
	U.resize(0);
	for(i=0;i<vertices;i++){
		U.push_back(0.0);
	};
	
	// desired angle sums are 2pi
	Theta.resize(0);
	for(i=0;i<vertices;i++){
		Theta.push_back(TWOPI);
	};

	relabel_vertices(infinity-1);		// move specified vertex to infinity
	
	compute_OPP();			// generate OPP data
	compute_TURN();			// generate TURN data
	test_consistency();		// check we have an honest triangulation with well-defined edge lengths
	cout << "File read and verified for consistency. " << vertices << " vertices.\n";
	cout << "time: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";
};


void packing::project_to_sphere(){		// project center list to sphere
	int i;
	center C;
	dbl X,Y;
	for(i=0;i<(int) POS.size();i++){
		X=POS[i].x;
		Y=POS[i].y;
		if(POS[i].x==0.0 && POS[i].y==0.0 && POS[i].t==1.0){
			// do nothing; this is the point at infinity
		} else {
			C.x=2.0*X/(X*X+Y*Y+1.0);
			C.y=2.0*Y/(X*X+Y*Y+1.0);
			C.t=(X*X+Y*Y-1.0)/(X*X+Y*Y+1.0);
			POS[i]=C;
		};
	};
};

void packing::restore_correct_labels(){		// restore labels from SWAP_LIST
	int i;
	center C;
	for(i=(int) SWAP_LIST.size();i>0;i=i-2){
		C=POS[SWAP_LIST[i-1]];
		POS[SWAP_LIST[i-1]]=POS[SWAP_LIST[i-2]];
		POS[SWAP_LIST[i-2]]=C;
	};
};

void packing::output_dual_format(ofstream &output_file){		// write list of centers to file
	int i;

	output_file << "[";
	for(i=0;i<(int) POS.size();i++){
		output_file << "[" << POS[i].x << "," << POS[i].y << "," << POS[i].t << "],\n";
	};
	output_file << "\"END\"];\n";
};

/*

The ideal output is
-----------------------------
[[x1,y1,z1],
...
[xn,yn,zn],
"END"];
-----------------------------
giving the positions on the sphere S^2 in R^3 of the points 1...n.

*/
