/* lobachevsky.cc definition of Lobachevsky function; 
	also some elementary vector and trig functions */
	
dvec operator-(dvec A, dvec B){
	dvec C;
	C.resize(0);
	int i;
	for(i=0;i<(int) A.size();i++){
		C.push_back(A[i]-B[i]);
	};
	return(C);
};

dbl norm(dvec A){
	dbl t;
	t=0.0;
	int i;
	for(i=0;i<(int) A.size();i++){
		t=t+(A[i]*A[i]);
	};
	return(t);
};

dbl dot(dvec A, dvec B){
	dbl t;
	t=0.0;
	int i;
	for(i=0;i<(int) B.size();i++){
		t=t+(A[i]*B[i]);
	};
	return(t);
};

dbl L_to_length(dbl L){	// length is e^{L/2}
	return(exp(L/2.0));
};

dbl angle(dbl r, dbl s, dbl t){	
	/* for Euclidean triangle with edge lengths r,s,t 
	returns angle opposite edge of length r	*/
	dbl alpha,u;
	
	// r^2 = s^2+t^2-2stcos(alpha)
	u=(s*s+t*t-r*r)/(2.0*s*t);
	if(u>=1.0){
		return(0.0);
	} else if (u<=-1.0){
		return(PI);
	} else {
		alpha = acos((s*s+t*t-r*r)/(2.0*s*t));	// warning: unstable near 1 or -1
		return(alpha);
	};
};

dbl Lambda(dbl x){
	// returns Lambda(x):= -\int_0^x \log{|2\sin(t)|} dt
	dbl t,l;
	dbl mesh_size;
	
	mesh_size=0.00001;
	l=0.0;
	for(t=mesh_size;t<x;t=t+mesh_size){
		l=l-(log(abs(2.0*sin(t)))*mesh_size);
	};
	
	return(l);
};

/* Note: computation of Lambda is reasonably accurate but slow. Instead, we should
  precompute Lambda at a fine mesh, and then when we need to evaluate it, linearly
  interpolate from precomputed values.	*/