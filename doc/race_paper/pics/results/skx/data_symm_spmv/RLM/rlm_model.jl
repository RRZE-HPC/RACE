#!/usr/bin/env julia

using StatsBase
const prgm = basename(Base.source_path())

if length(ARGS) < 5
   println("Usage: ", prgm, " bw table_nnz table_spmv table_measured_alpha filename")
   exit(1)
end

bw=parse(Float64,ARGS[1])
table_nnz=ARGS[2]
table_spmv=ARGS[3]
table_measured_alpha=ARGS[4]
filename=ARGS[5]

nnz=readdlm(table_nnz)

names=nnz[1:end,1]
nnzr=nnz[1:end,2]

spmv=readdlm(table_spmv)

spmv_names=spmv[1:end, 1]
spmv_perf=spmv[1:end, 2]

measured_alpha=readdlm(table_measured_alpha)

measured_alpha_names=measured_alpha[1:end,1]
measured_alpha_val=measured_alpha[1:end,2]

max_width=0
for i in 1:length(names)
	max_width=max(max_width, strwidth(names[i]))
end
#Now create new string

for i in 1:length(names)
	curr_width=strwidth(names[i])	
	pad_str=""
	for j in curr_width:max_width+1
		pad_str=string(pad_str, " ")
	end
	names[i]=string(names[i], pad_str)
end

		
function calc_alpha_from_spmv(id, bw)
	perf=spmv_perf[id]
	return ((2*bw/perf)-(12+20/nnzr[id]))/8.0
end

function perf_symm_spmv_from_spmv_perf(id, bw, opt)
	nnzr_symm=(nnzr[id]-1)/2+1
	alpha = 0
	if(opt)
		alpha = 1/nnzr_symm
	else 
		alpha=calc_alpha_from_spmv(id, bw)
		if(alpha < (1/nnzr_symm))
			info("For matrix ", spmv_names[id], " alpha = ", alpha, " but opt alpha =", 1/nnzr_symm)
			alpha = 1/nnzr_symm #reset it
		end
	end
	return (4/(8+4+24*alpha+(4)/nnzr_symm))*bw
end

alphas=zeros(length(nnzr))
for i in 1:length(nnzr)
	nnzr_symm=(nnzr[i]-1)/2+1
	alpha=calc_alpha_from_spmv(i, bw)
	if(alpha < (1/nnzr_symm))
		info("For matrix ", spmv_names[i], " alpha = ", alpha, " but opt alpha =", 1/nnzr_symm)
		alpha = 1/nnzr_symm #reset it
	end
	alphas[i]=alpha
end

measured_alphas=zeros(length(nnzr))
for i in 1:length(nnzr)
	measured_alphas[i]=measured_alpha_val[i]
end

function perf_symm_spmv_from_measured_alpha(id, bw)
	nnzr_symm=(nnzr[id]-1)/2+1
	alpha = 0
	alpha = measured_alpha_val[id]
	return (4/(8+4+24*alpha+(4)/nnzr_symm))*bw
end

id=Array{Int64,1}(length(nnzr))
for i in 1:length(nnzr)
	id[i] = i
end

symm_spmv_perf_alpha=perf_symm_spmv_from_spmv_perf.(id,bw,false)
symm_spmv_perf_opt=perf_symm_spmv_from_spmv_perf.(id,bw,true)
symm_spmv_measured_alpha=perf_symm_spmv_from_measured_alpha.(id,bw)

deviation=100*(alphas-measured_alphas)./measured_alphas

writedlm(string(filename),[0:length(names)-1 names nnzr symm_spmv_measured_alpha symm_spmv_perf_opt symm_spmv_perf_alpha],"|")
writedlm(string("alpha.txt"),[0:length(names)-1 names nnzr measured_alphas alphas deviation],",\t")
