using DifferentialEquations, Random, MAT, NaNMath, Statistics
include("integrator_ern.jl")

# 1 is restarting
# use this during self assembly to create burster, then turn off homeostatic mechanism
# when homeostatic mechanism if off, this burster may change
restarter = 1;

tgt1 = 0.05 
tgt2 = 0.07 
tgt3 = 0.3 
gf = 10.0
gs = 3.0

function extract_last_values(filepath::String)
    last_values = []

    for line in eachline(filepath)
        if occursin("Last Value", line)
            # Extract the part after the '|'
            parts = split(line, "|")
            if length(parts) >= 2
                array_str = strip(parts[end])
                # Remove square brackets and split into numbers
                array_clean = strip(array_str, ['[', ']'])
                values = parse.(Float64, split(array_clean))
                push!(last_values, values)
            end
        end
    end

    return last_values
end
data = extract_last_values("/Users/linxiao/Desktop/july.29 simulation/random_ini_0.07_0.1_0.5_2.txt")
data_matrix = hcat(data)
# Optionally convert to matrix if all vectors are the same length
print(data_matrix)
# for i= 1:1:1
for (i, IC) in enumerate(data)
	ICs=reshape(IC, 1, :)
	
	# chunkSize = 3000.0;  # this is in seconds! (not milliseconds)
	# tEnd = chunkSize*29; # modify the total time here
	
	num=i
	hsh = randstring(['A':'Z'; '0':'9'], 6);
	# path = "/Users/linxiao/Desktop/July.24 Simulation/final bursting_txt7/"
	path = "/Users/linxiao/Desktop/july.29 simulation/0.07_0.1_0.5_finalbursting_2/"
	# fname = path*"tgt1."*string(tgt1)*".tgt2."*string(tgt2)*".tgt3."*string(tgt3)*"."*hsh;
	# mkpath(fname)

	# Initial condition
	# ICs = [-50 0.004 0.26 0.02 0.357 0.97 0.71 25.09 11.2 0.03 0.42 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0]; 
	# ICs = [-50 0.004 0.26 0.02 0.357 0.97 0.45 20 9 0.05 0.45 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];
	
	#print("\n RUN NUMBER: ",i, "\n")
	# alphaVal = integrator_v1(fname,tgt1,tgt2,tgt3,gf,gs,ICs,tEnd,0.0001,chunkSize,[300000.0, 300000.0, 150000.0])

	
	## restart to capture burster w/o homeostasis
	if restarter == 1
		chunkSize = 60; tEnd = chunkSize*1; # modify the total time here
		## LOAD ICs previously run simulations
		# file = matopen(fname*"/"*hsh*"_029.mat"); # if you modify the total time, remember to change the name of the file.
		# ICs = read(file,"lastValue")
		# close(file)
		##
		## rename for new run
		fname = path*"final"*string(i)*".RESTART."
		mkdir(fname)
		#1
		# ICs = [-55.8073084228253 0.0012271214028683978 0.5527154630024492 0.018170056976805196 0.763661879712376 0.8755098736756551 0.421944294055214 39.18058913292683 18.674924089392924 0.060277953670564904 3.4833853384119635 0.04033205016635652 0.9019873067357356 0.15821753094989618 0.5764683978156114 0.5378060698174166 0.023330591563658597 -0.0033718788997923192 0.020608613697997217 0.000966789014908481]
		#2
		# ICs=[-54.78985178969181 0.0015819624246156555 0.5398028515343072 0.021114380047614175 0.9691926857106165 0.916985887351939 0.5350244767362176 5.851435457678389 1.2323567275874783 0.0396316459285954 3.8054820478785785 0.0334370331192717 0.9178982445126929 0.1339899420023511 0.6234231441243723 0.5266445344551137 -0.02441737430124225 0.013982857738835458 -0.01197519623507591 0.00015667226908518624]

		alphaVal = integrator_v1(fname,tgt1,tgt2,tgt3,gf,gs,ICs,tEnd,0.0001,chunkSize,[Inf, 300000.0, 150000.0])
	end 
end