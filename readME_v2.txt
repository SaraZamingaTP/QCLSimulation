2023.04.17 - Sara's read me file

Updated project: this folder contains two versions for each .m file:
1) First version: FP laser;
2) Second version: DFB laser;

	a) DeviceData_OptExpr_DFB: file for parameters setting; with respect to the FP
		version, there is the possibility to externally set the alpha factor and kDFB. 
	b) Main_TDTW_DFBQCL: file to run simulations: it performs alpha factor and k sweeps,
		according to the input flags isAlphaVar, isKcoupVar; 
		if both false, it just runs the simulation for the input I current value
	c) TDTW_runge_DFBQCL: file for algorithm execution: solution of field's equations,
		boundary conditions, carriers' equation. 
	d) post_proc_simple_NEW: - updated figure settings;
				 - representation of the optical spectrum (dB) as a function
				of lambda (not frequency, GHz)
	e) ParallelSim: file I use to perform multiple simulations simultaneously;
			isAlphaVar, isKcoupVar are two flags used to determine when to
			execute the parfor loop. 
	f) RemoteParallelSim: file to perform the simulations using the RemoteServer. 
	g) LIcurve: function used to represent P(I)
	h) OpticalSpectrum: function to represent the optical spectrum as a function of
				lambda
	i) PowerSpectrumMap: function to represent the optical spectrum as a function 
				of the current

-- Results: folder containing simulations' results