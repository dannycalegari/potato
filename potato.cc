/************************************
*                                   *
*	potato.cc version 0.01			*
*	September 17 2012				*
*									*
*									*
*	Copyright Danny Calegari 2012	*
*	Released under the GPL license	*
*									*
************************************/

// standard libraries to include	

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <complex>
#include <sstream>
#include <ctime>
#include <assert.h>
#include <limits>



using namespace std;

// preprocessor definitions

#define PI 			3.141592653589793238462643383279
#define TWOPI 		6.283185307179586476925286766558
#define PIONTWO		1.5707963267948966192313216916395
#define ivec		vector<int>
#define imat		vector<vector<int> >
#define ivl			vector<vector<int> >	// list of ivecs
#define dbl			long double					// can replace this with arbitrary precision later
#define dvec		vector<dbl>
#define dmat		vector<vector<dbl> >
#define dvl			vector<vector<dbl> >	// list of dvecs

#include "graphics.cc"
#include "lobachevsky.cc"
#include "vector.cc"
#include "packing.cc"
#include "input_output.cc"
#include "dual_format.cc"
#include "layout.cc"

// global constants

typedef std::numeric_limits< dbl > dubbl;


int main(int argc, char *argv[]){
	packing P;
	ifstream input_file;
	ofstream output_file;
	string S;
	bool b,finished;
	
	cout.precision(dubbl::digits10);

	if(argc>1){
		S=argv[1];
		if(S=="-nl"){
			input_file.open(argv[2]);
			P.input_data(input_file,false);
		} else if(S=="-df"){
			input_file.open(argv[2]);
			P.dual_format_input_data(input_file);
		} else {
			input_file.open(S.c_str());
			P.input_data(input_file,true);
		};
		input_file.close();
	} else {
		cout << "Enter name of input file: ";
		cin >> S;
		cout << "Variable lengths? (1 for yes, 0 for no): ";
		cin >> b;
		input_file.open(S.c_str());	// function .c_str() converts string to char*
		P.input_data(input_file,b);
		input_file.close();
	};	
		
	P.flow_to_packing(0.0000000001);
	P.test_consistency();
	P.determine_layout_order();
	P.determine_layout_position();
	
	setup_graphics();
	usleep(100000);
	P.draw_packing();
	P.draw_eps();
	XFlush(display);
	finished=false;
	while(finished==false){ 
		XNextEvent(display, &report);
        switch (report.type) {
        	case KeyPress:
				finished=true;
				XCloseDisplay(display);
				exit(0);
				break;
			default:
				break;
		};
	};
	
	if(S=="-df"){
		P.project_to_sphere();
		P.restore_correct_labels();
		output_file.open("center_list.txt");
		output_file.precision(dubbl::digits10);
		P.output_dual_format(output_file);
		output_file.close();
		return 0;
	};

	
	return 0;
}
