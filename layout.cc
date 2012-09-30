/* layout.cc  layout functions */

/*  Layout.  By convention, vertex 0 is the "center", and ATT[0]=0. Then for all other i
	except the "infinite vertex", ATT[i] is some j which is closer to vertex 0 than i.
	Hence i -> ATT[i] determines a deformation retraction to vertex 0. 
	We also insist that ATT[i] is not a boundary vertex.
	
	By convention, LAY[0]=0 and for all other i, LAY[i] > LAY[ATT[i]].

		ivec LAY;	// layout order
		ivec ATT;	// attached to
		vector<center> POS;		// layout position	*/
		
void packing::determine_layout_order(){		// compute ATT and LAY
	int i,j,k,v;
	bool all_found;
	clock_t start;
	start=clock();
	
	cout << "determining ATT and LAY vectors.\n";

	v=ADJ[(int) ADJ.size()-1].size()+1;	// v is number of boundary vertices + 1
	cout << "internal vertices are 0 to " << ADJ.size()-1-v << "\n";
	cout << "boundary vertices are " << ADJ.size()-v << " to " << ADJ.size()-2 << "\n";
	cout << "infinite vertex is " << ADJ.size()-1 << "\n";
	all_found=false;	// initial value
	ATT.resize(0);
	LAY.resize(0);
	for(i=0;i<(int) ADJ.size();i++){
		ATT.push_back(-1);
	};
	ATT[0]=0;
	LAY.push_back(0);
	while(all_found==false){
		k=0;
		for(i=0;i<(int) ADJ.size()-1;i++){	// for each vertex i *except the infinite vertex*
			if(ATT[i]==-1){		// do we already know ATT[i]?
				k++;		// if not, notice this
				for(j=0;j<(int) ADJ[i].size();j++){		// for each neighbor of i
					if(ADJ[i][j]<=(int) ADJ.size()-1-v && ATT[ADJ[i][j]]!=-1){		
					// if we know ATT for the neighbor AND the neighbor is not a boundary vertex
						ATT[i]=ADJ[i][j];		// set ATT[i] to the neighbor.
						LAY.push_back(i);		// add to layout queue
						j=(int) ADJ[i].size();	// skip to the end
					};
				};
			};
		};
		if(k==0){
			all_found=true;
		};
		cout << k << " left to determine.\n";
	};
	cout << "time: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";
};

void packing::determine_layout_position(){		// compute POS
	int i,j,k,l,m;
	center C,D;
	dbl a,b;
	clock_t start;
	start=clock();
	
	cout << "Determining positions. \n";
	
	C.x=0.0;
	C.y=0.0;
	C.t=0.0;
	
	POS.resize(0);
	for(i=0;i<(int) ADJ.size();i++){
		POS.push_back(C);		// POS[i]=(0,0,0) is default
	};
	POS[ADJ.size()-1].t=1.0;

	for(i=1;i<(int) LAY.size();i++){
		j=ATT[LAY[i]];	// LAY[i] is the ith to be laid out, and j is the vertex it rests on
		D=POS[j];
		for(k=0;k<(int) ADJ[j].size();k++){
			if(ADJ[j][k]==LAY[i]){
				l=k;
			};
		};
		// at this point ADJ[j][l]=LAY[i] and we know POS[j].
		assert(ADJ[j][l]==LAY[i]);
		a=D.t;	// a is angle of edge j -0->
		for(k=0;k<l;k++){
			a=a+Atilde(j,k);
		};
		// at this point a is angle of edge j -l-> i
		C.x=D.x+cos(a)*L_to_length(Ltilde(j,l));
		C.y=D.y+sin(a)*L_to_length(Ltilde(j,l));
		// C.x,C.y is location of vertex i.
		m=OPP[j][l];
		assert(ADJ[LAY[i]][m]==j);	// check that ADJ[LAY[i]][m]=j
		b=0.0;
		for(k=0;k<m;k++){
			b=b+Atilde(LAY[i],k);
		};
		// at this point b+C.t=a+PI so C.t = a-b+PI
		C.t=a-b+PI;
		while(C.t<0.0){
			C.t=C.t+TWOPI;
		};
		while(C.t>TWOPI){
			C.t=C.t-TWOPI;
		};
		POS[LAY[i]]=C;	// set position LAY[i].
	};
	rescale();
	cout << "time: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "\n";
};

void packing::rescale(){		// rescale so radius of 0 and radius of infinity are inverse
	dbl outer,t;
	int i,j;
	
	outer=0.0;
	for(j=0;j<(int) ADJ[(int) ADJ.size()-1].size();j++){
		i=ADJ[(int) ADJ.size()-1][j];
		t=POS[i].x*POS[i].x + POS[i].y*POS[i].y;
		if(outer<t){
			outer=t;
		};
	};
	// rescale by 1.0/outer*inner
	t=1.0/(outer);
	t=sqrt(t);
	for(i=0;i<(int) ADJ.size();i++){
		POS[i].x=POS[i].x*t;
		POS[i].y=POS[i].y*t;
	};
};


void packing::draw_packing(){		// graphical output
	int i,j;
	point P,Q;
	for(i=0;i<(int) ADJ.size()-1;i++){
		P=center_to_point(POS[i]);
	//	draw_circle(P,3,0x000000);
		for(j=0;j<(int) ADJ[i].size();j++){
			if(i<ADJ[i][j] && (ADJ[i][j]!=(int) ADJ.size()-1)){
				Q=center_to_point(POS[ADJ[i][j]]);
				draw_line(P,Q,0x000000);		
			};
		};
	};
	XFlush(display);
};

void packing::draw_eps(){		// eps output
	int i,j;
	ofstream eps_out_file;
	
	typedef std::numeric_limits< dbl > dubbl;

	
	eps_out_file.open("potato_packing.eps");
	eps_out_file.precision(dubbl::digits10);

    eps_out_file << "%!PS-Adobe-2.0 EPSF-2.0 \n";
    eps_out_file << "%%BoundingBox: 0 0 500 500 \n";
    eps_out_file << "gsave 200 200 scale 1.25 1.25 translate 0 setlinewidth \n";
    eps_out_file << "/l {4 dict begin /y2 exch def /x2 exch def /y1 exch def /x1 exch def \n";
    eps_out_file << "newpath x1 y1 moveto x2 y2 lineto stroke end} def \n";
//	point P,Q;
	for(i=0;i<(int) ADJ.size()-1;i++){
//		P=center_to_point(POS[i]);
		for(j=0;j<(int) ADJ[i].size();j++){
			if(i<ADJ[i][j] && (ADJ[i][j]!=(int) ADJ.size()-1)){
//				Q=center_to_point(POS[ADJ[i][j]]);
				eps_out_file << POS[i].x << " " << POS[i].y << " " << POS[ADJ[i][j]].x << " " << POS[ADJ[i][j]].y << " l\n";
			};
		};
	};
	eps_out_file << "grestore %eof \n";
};

