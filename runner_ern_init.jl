using DifferentialEquations, Random, MAT, NaNMath, Statistics
include("integrator_ern_test.jl")
# resultsFilePath="/Users/linxiao/Desktop/July.21 simulation/convergence_results2.txt"
tgt1 = 0.07
tgt2 = 0.1
tgt3=0.5
gf = 10.0
gs = 3.0
fname="/Users/linxiao/Desktop/July.29 simulation/random_ini_0.07_0.1_0.5_2.txt"
global convergence_count=0
convergence_set_initial=[]
convergence_set_final=[]
# mkpath(fname)
for i=1:1:50
	chunkSize = 3000.0;  #this is in seconds! (not milliseconds) 
	tEnd = chunkSize*29; # total time is 29 * 3000 seconds
	hsh = randstring(['A':'Z'; '0':'9'], 6);
	# gNap = round(rand() * (0.99 - 0.16) + 0.16, digits=2)
	# gNa  = round(rand() * (29.75 - 20.49) + 20.49, digits=2)
	# gK   = round(rand() * (13.91 - 9.01) + 9.01, digits=2)
	# gCa  = round(rand() * (0.10 - 0.01) + 0.01, digits=2)
	# gCan = round(rand() * (1.00 - 0.01) + 0.01, digits=2)
	
	#basline 
	# gNap = round(rand(Float64), digits=2)  
	# gNa  = round(rand(Float64)*10 + 20, digits=2) 
	# gK   = round(rand(Float64)*5 + 9, digits=2)
	# gCa  = round(rand(Float64)*0.1, digits=2)
	# gCan = round(rand(Float64), digits=2)  

	#restricted by 32/50 DB
	gNap = round(rand() * (0.99 - 0.05) + 0.05, digits=2)
	gNa  = round(rand() * (29.76 - 20.28) + 20.28, digits=2)
	gK   = round(rand() * (13.61 - 9.11) + 9.11, digits=2)
	gCa  = round(rand() * (0.1 - 0.01) + 0.01, digits=2)
	gCan = round(rand() * (1.0 - 0.06) + 0.06, digits=2)

	#restricted by 40/50 range PB
	# gNap = round(rand(Float64)*(1.0 - 0.02) + 0.02, digits=2)
	# gNa  = round(rand(Float64)*(29.62 - 20.01) + 20.01, digits=2)
	# gK   = round(rand(Float64)*(13.99 - 9.13) + 9.13, digits=2)
	# gCa  = round(rand(Float64)*(0.1 - 0.01) + 0.01, digits=2)
	# gCan = round(rand(Float64)*(0.99 - 0.02) + 0.02, digits=2)



				# path = "/Users/linxiao/Desktop/July.21 simulation/S_300_alpha_150/"
				# fname = path*"tgt1."*string(tgt1)*".tgt2."*string(tgt2)*".tgt3."*string(tgt3)*"."*hsh;

	# ICs at the resting state, no opioid
	ICs = [-50 0.004 0.26 0.02 0.357 0.97 gNap gNa gK gCa gCan 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];

	# ICs=[-4.94602754e+01 5.95188454e-03 3.48892895e-01 3.68049644e-02 1.00037129e+00 8.67459910e-01 5.24467663e-01 1.63382865e+02 5.45666244e+01 1.04894125e-01 4.37130731e+00 3.96818078e-01 3.71084512e-01 6.95229200e-01 1.50931533e-01 8.13214094e-01 2.67230728e-02 7.94075960e-03 -1.90665175e-02 1.06169181e-01
	# ]

	# Starting at a particular point where the homeostatic mechanism turns off, opioids on 
	# ICs = [-5.71033206e+01  8.87827789e-04 8.07960285e-01 2.49639611e-02 9.62987480e-01  9.39178059e-01 4.41431632e-01  1.43670990e+02 3.34950769e+01 8.82863765e-02 1.35093399e+00 4.73948267e-02 8.86021351e-01 1.82210745e-01 5.32849206e-01 4.98078781e-01 3.30733918e-01 9.84100317e-02 1.57550906e-01 8.64526022e-04 ]
	
	#print("\n RUN NUMBER: ",i, "\n")
	alphaVal = integrator_v1(fname,convergence_count,convergence_set_initial,convergence_set_final,tgt1,tgt2,tgt3,gf,gs,ICs,tEnd,0.0001,chunkSize,[300000.0, 300000.0, 150000.0])
	end
open(fname, "a") do io
	println(io, "convergence count=$convergence_count")
	println(io, "convergence_set_initial=$convergence_set_initial")
	println(io, "convergence_set_final=$convergence_set_final")
end		