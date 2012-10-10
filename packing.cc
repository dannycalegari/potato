/* packing.cc	*/

struct center{
	/* specifies location of a vertex. x and y are real and imaginary parts, t
	is the angle that edge 0 is rotated from the positive real axis. 
	Later, when we project stereographically, x,y,t are x,y,z coordinates.	*/

	dbl x,y,t;
};

point center_to_point(center C){
	point P;
	P.x = 400+(int) 300.0*C.x;
	P.y = 400-(int) 300.0*C.y;
	return(P);
};

class packing{
	public:
		// data defining discrete conformal metric

/*	Notation for vertex is an integer i. Notation for edge emanating from i is i -j->. 		
	ADJ[i][j] is the jth neighbor of vertex i, so if ADJ[i][j]=k we have i -j-> k.
	If ADJ[i][j]=k and OPP[i][j]=l then k -l-> i; i.e. ADJ[j][OPP[i][j]]=i.
	If ADJ[i][j]=k and TURN[i][j]=l then k -l-> m; if TURN[k][l]=n then m -n-> i.
	L[i][j] is 2 log length(i -j-> k) */
	
		ivl	ADJ;	// adjacency list
		dvl	L;		// edge log lengths
		ivl OPP;	// data to recover reverse of edge
		ivl TURN;	// data to recover next edge of triangle
		
		// functions to precompute after reading in adjacency data
		void compute_OPP();		// compute OPP from ADJ
		void compute_TURN();	// compute TURN from ADJ
		void test_consistency();	// test data for consistency
		ivec SWAP_LIST;		// list of swaps performed on vertices, in order. Undo to make labels consistent at the end.
		void swap_vertices(int,int);	// swaps labels of two vertices
		void relabel_vertices(int);		// sets specific vertices and its neighbors to the end
		
		// input/output functions
		
/* functions read/write adjacency list and edge log lengths from/to specified files */
		
		void input_data(ifstream &, bool);	// read from file in standard format
			// if bool==true log lengths of edges are specified; otherwise they are all set to 0.
		void dual_format_input_data(ifstream &);		// read from file in dual format
		void output_data(ofstream &);	// write to file
		void write();					// write to cout
		
/*	U is the vector of weights at the vertices; U[i] is the weight at vertex i.	
	Theta is the desired angle sums at vertices; Theta[i] is desired sum at vertex i. */
		
		// variables relevant to packing
		dvec U;			// vertex weights
		dvec Theta;		// desired angle sums at vertices
		
/*  edge and angle functions; implicitly, functions of U.	*/

		dbl Ltilde(int, int);	// adjusted edge log length
		dbl Atilde(int, int);	// angle at vertex in triangle
		dbl E();				// functional E
		dbl Fitness(int);		// how far do angle sums deviate from Theta?
		dvec JAC();				// dE/dU
		dmat DJAC();			// d^2E/dU^2
		
/*  flow to packing	*/

		void flow_to_packing(dbl);	// adjust U until angle sums at all interior vertices equal Theta with prescribed error
		void Newton_packing(dbl);	// adjust U by Newton's method

/*  Layout.  By convention, vertex 0 is the "center", and ATT[0]=-1. Then for all other i
	except the "infinite vertex", ATT[i] is some j which is closer to vertex 0 than i.
	Hence i -> ATT[i] determines a deformation retraction to vertex 0. 
	
	By convention, LAY[0]=0 and for all other i, LAY[i] > LAY[ATT[i]]. */

		ivec LAY;	// layout order
		ivec ATT;	// attached to
		vector<center> POS;		// layout position
		void determine_layout_order();		// compute ATT and LAY
		void determine_layout_position();	// compute POS
		void draw_packing();	// graphical output
		void draw_eps();	// eps output file
		void rescale();		// rescale so radius of 0 and radius of infinity are inverse
		void project_to_sphere();	// project center list to sphere
		void restore_correct_labels();	// restore labels from SWAP_LIST
		void output_dual_format(ofstream &);	// write list of centers to file
};

void packing::compute_OPP(){
	OPP.resize(0);
	ivec O;
	int i,j,k,l;
	for(i=0;i<(int) ADJ.size();i++){
		O.resize(0);
		for(j=0;j<(int) ADJ[i].size();j++){
			k=ADJ[i][j];
			for(l=0;l<(int) ADJ[k].size();l++){
				if(ADJ[k][l]==i){
					O.push_back(l);
					l=(int) ADJ[k].size();
				};
			};
		};
		OPP.push_back(O);
	};
};

void packing::compute_TURN(){
	TURN.resize(0);
	ivec T;
	int i,j,k,l;
	for(i=0;i<(int) ADJ.size();i++){
		T.resize(0);
		for(j=0;j<(int) ADJ[i].size();j++){
			k=ADJ[i][j];
			l=OPP[i][j];
			if(l==0){
				l=(int) ADJ[k].size()-1;
			} else {
				l=l-1;
			};		
			T.push_back(l);
		};
		TURN.push_back(T);
	};
};

void packing::test_consistency(){
	int i,j,k,l,m,n;
	for(i=0;i<(int) ADJ.size();i++){
		for(j=0;j<(int) ADJ[i].size();j++){
			k=ADJ[i][j];
			l=OPP[i][j];
			assert(ADJ[k][l]==i);
			assert(L[k][l]==L[i][j]);
			assert(Ltilde(k,l)==Ltilde(i,j));	// order of addition matters for roundoff
			l=TURN[i][j];
			m=ADJ[k][l];
			n=TURN[k][l];
			assert(ADJ[m][n]==i);
		};
	};
};

void packing::swap_vertices(int k, int l){	// swaps labels of two vertices
	int i,j;
	ivec S;
	dvec T;
	dbl u;
	// switch pointers to k and to l from other vertices
	for(i=0;i<(int) ADJ.size();i++){
		for(j=0;j<(int) ADJ[i].size();j++){
			if(ADJ[i][j]==k){
				ADJ[i][j]=l;
			} else if(ADJ[i][j]==l){
				ADJ[i][j]=k;
			};
		};
	};
	S=ADJ[k];		// swap labels at k and l
	ADJ[k]=ADJ[l];
	ADJ[l]=S;
	
	T=L[k];			// swap lengths at k and l
	L[k]=L[l];
	L[l]=T;
	
	u=U[k];			// swap vertex weights
	U[k]=U[l];
	U[l]=u;
	SWAP_LIST.push_back(k);
	SWAP_LIST.push_back(l);
};


void packing::relabel_vertices(int k){
	/* sets vertex k and its neighbors to the ``last'' vertices
	Effectively sends k ``to infinity'' and its neighbors to boundary vertices	*/
	int i,l,m;
	l=(int) ADJ.size()-1;
	swap_vertices(k,l);	// swap k for last vertex
	m=ADJ[l].size();
	for(i=0;i<(int) m;i++){	// swap neighbors of last for last-m . . last-1
		swap_vertices(ADJ[l][i],l-m+i);
	};
	compute_OPP();
	compute_TURN();
	test_consistency();
};

dbl packing::Ltilde(int i, int j){	// adjusted log edge length
	int k;
	k=ADJ[i][j];
	if(i<k){
		return(L[i][j]+U[i]+U[k]);
	} else {
		return(L[i][j]+U[k]+U[i]);
	};
};

dbl	packing::Atilde(int i, int j){	// angle at i between edges j and j+1 in adjusted metric
	dbl r,s,t;	// edges of triangle are s,r,t; we want angle opposite r.
	dbl alpha;
	int k,l,m,n;
	k=ADJ[i][j];
	l=TURN[i][j];
	m=ADJ[k][l];
	n=TURN[k][l];
	// s,r,t are adjusted lengths
	s=L_to_length(Ltilde(i,j));
	r=L_to_length(Ltilde(k,l));
	t=L_to_length(Ltilde(m,n));
	alpha=angle(r,s,t);
	
	return(alpha);
};

dbl packing::Fitness(int k){
	/* sum of absolute value of angle defects at interior vertices
	   note: we assume that the non-interior vertices are the last k values!
	   this can be achieved by first doing packing::relabel_vertices(i) 
	   where i is the vertex to move ``to infinity'' */
	   
	int i,j;
	dbl e,A;
	e=0.0;
#pragma omp parallel for private(A) reduction(+:e)
	for(i=0;i<(int) ADJ.size()-k;i++){	// 7
		A=0.0;
		for(j=0;j<(int) ADJ[i].size();j++){
			A=A+Atilde(i,j);
		};
		e=e+abs(A-Theta[i]);	// angle deviation
	};
	return(e);
};

dbl	packing::E(){	// functional
	dbl e;
	int i,j;
	e=0.0;
#pragma omp parallel for private(j) reduction(+:e)
	for(i=0;i<(int) ADJ.size();i++){
		for(j=0;j<(int) ADJ[i].size();j++){
			e=e+Ltilde(i,j)*(Atilde(i,j)-PIONTWO);
			e=e+2.0*Lambda(Atilde(i,j));	// costly to evaluate right now
		};
		e=e+U[i]*Theta[i];
	};
	
	return(e);
};

dvec packing::JAC(){	// dE/dU
	// Note: ignores boundary and infinite vertices
	dvec E;
	dbl e;
	int i,j;
	int v,vv;
	
	v=(int) ADJ.size()-1;
	vv=1+ADJ[v].size();
	
	E.resize(0);
	dbl Edata[ADJ.size()];
#pragma omp parallel for private(j,e)
	for(i=0;i<(int) ADJ.size()-vv;i++){
		e=Theta[i];
		for(j=0;j<(int) ADJ[i].size();j++){
			e=e-Atilde(i,j);
		};
//		E.push_back(e);
		Edata[i] = e;
	};
	for(i=0;i<(int) ADJ.size()-vv;i++)
		E.push_back(Edata[i]);
	return(E);
};	

dmat packing::DJAC(){	// d^2E/dU^2
	// Note: ignores boundary and infinite vertices
	dmat M;
	dbl e;
	dvec E;
	int i,j,k,l;
	int v,vv;
	
	v=(int) ADJ.size()-1;
	vv=1+ADJ[v].size();
	
	M.resize(0);
	for(i=0;i<(int) ADJ.size()-vv;i++){
		E.resize(0);
		for(j=0;j<(int) ADJ.size()-vv;j++){
			for(k=0;k<(int) ADJ[i].size();k++){	// is j a neighbor of i?
				e=0.0;	// default value
				if(ADJ[i][k]==j){	 // if yes,
					l=OPP[i][k];	// ADJ[j][l]=i
					e=-(1.0/tan(Atilde(i,k)))-(1.0/tan(Atilde(j,l)));
				};
			};
			E.push_back(e);
		};
		M.push_back(E);
	};
	return(M);
};

void packing::flow_to_packing(dbl accuracy){	// adjust U until angle sums at all vertices equal Theta
	dbl speed, fit, ang, bound;
	dvec J,JOLD;
	dvec S,SOLD;
	int i,j,v,count;
	clock_t start;
	start=clock();
	cout << "Solving packing problem.\n";
	/* Vertex v at infinity has j-1 boundary vertices.
	Consequently we should not adjust U for the last j vertices, or count the
	angle defects there in the contribution to Fitness. */
	v=(int) ADJ.size()-1;
	j=1+ADJ[v].size();
	
	count=0;
	fit=Fitness(j);
	bound=10.0;

	while(fit>accuracy){			// sum of abs of angle defects for all but last j vertices
//		fit=Fitness(j);
		J=JAC();
		S=U;
		if(count==0){		// should determine first step size by line search
			speed=0.001;
		} else {
			speed=dot(S-SOLD,J-JOLD)/norm(J-JOLD);	// Barzilai-Borwein gradient method
		};
		fit=0.0;
#pragma omp parallel for private(i)		// let's see what this does!
		for(i=0;i<(int) U.size()-j;i++){	// adjust U for all but last j vertices
			ang=J[i];
			fit=fit+abs(ang);
			U[i]=U[i]-ang*speed;
		};
		JOLD=J;
		SOLD=S;
		if(fit<bound){
//		if(count%100==0){
			cout << "After " << count << " steps fitness is " << fit << "\n";
			bound=bound*0.25;
		};
		count++;
	};
	cout << "time: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";
};
