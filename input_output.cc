/* input_output.cc	input output routines */

void write(dvec J){
	int i;
	for(i=0;i<(int) J.size();i++){
		cout << J[i] << "\n";
	};
	cout << "\n";
};

void write(ivec J){
	int i;
	for(i=0;i<(int) J.size();i++){
		cout << J[i] << " ";
	};
	cout << "\n";
};

void packing::input_data(ifstream &input_file, bool variable_lengths){
	ivec I;	
	dvec T;	
	
	ADJ.resize(0);	// initialize data
	L.resize(0);
	SWAP_LIST.resize(0);
	
	int i,j,vertices,valence;
	int k;
	dbl l;
	
	input_file >> vertices;
	
	for(i=0;i<vertices;i++){
		I.resize(0);
		T.resize(0);
		input_file >> valence;
		for(j=0;j<valence;j++){
			input_file >> k;	// adjacent vertex
			if(variable_lengths){
				input_file >> l;	// log length of edge
			} else {
				l=0.0;
			};
			I.push_back(k);
			T.push_back(l);
		};
		ADJ.push_back(I);
		L.push_back(T);
	};
		
	// initial vertex weights are 0
	U.resize(0);
	for(i=0;i<vertices;i++){
		U.push_back(0.0);
	};
	
	// desired angle sums are 2pi
	Theta.resize(0);
	for(i=0;i<vertices;i++){
		Theta.push_back(TWOPI);
	};

//	relabel_vertices(vertices-1);	// adjust labels on neighbors of last vertex
	relabel_vertices(0);	// adjust labels on neighbors of last vertex
	
	compute_OPP();			// generate OPP data
	compute_TURN();			// generate TURN data
	test_consistency();		// check we have an honest triangulation with well-defined edge lengths
	cout << "File read and verified for consistency. " << vertices << " vertices.\n";
};

void packing::output_data(ofstream &output_file){
	int i,j;
	output_file << (int) ADJ.size() << "\n";
	for(i=0;i<(int) ADJ.size();i++){
		output_file << (int) ADJ[i].size() << "    ";
		for(j=0;j<(int) ADJ[i].size();j++){
			output_file << ADJ[i][j] << " ";
			output_file << L[i][j] << " ";
		};
		output_file << "\n";
	};
};

void packing::write(){
	int i,j;
	cout << "Data for packing.\n";
	cout << "Number of vertices is " << (int) ADJ.size() << "\n";
	for(i=0;i<(int) ADJ.size();i++){
		cout << "Valence of vertex " << i << " is " << (int) ADJ[i].size() << "\n";
		for(j=0;j<(int) ADJ[i].size();j++){
			cout << "Neighbor " << j << " is " << ADJ[i][j] << "  ";
			cout << "Length of edge " << i << "-" << ADJ[i][j] << " is " << L[i][j] << "\n";
			cout << "Adjusted length of edge " << i << "-" << ADJ[i][j] << " is " << Ltilde(i,j) << "\n";
		};
	};
	cout << "E is " << E() << "\n";
};