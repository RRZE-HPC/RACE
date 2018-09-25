#!/usr/bin/env julia

using StatsBase
const prgm = basename(Base.source_path())

if length(ARGS) < 4
   println("Usage: ", prgm, " bw table_nnz table_spmv filename")
   exit(1)
end

bw=parse(Float64,ARGS[1])
table_nnz=ARGS[2]
table_spmv=ARGS[3]
filename=ARGS[4]

nnz=readdlm(table_nnz)

names=nnz[1:end,1]
nnzr=nnz[1:end,2]

spmv=readdlm(table_spmv)

spmv_names=spmv[1:end, 1]
spmv_perf=spmv[1:end, 2]

function calc_alpha_from_spmv(id, bw)
	perf=spmv_perf[id]
	return ((2*bw/perf)-(12+16/nnzr[id]))/8.0
end


function perf_symm_spmv(id, bw, opt)
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

id=Array{Int64,1}(length(nnzr))
for i in 1:length(nnzr)
	id[i] = i
end

symm_spmv_perf_alpha=perf_symm_spmv.(id,bw,false)
symm_spmv_perf_opt=perf_symm_spmv.(id,bw,true)

writedlm(string(filename),[0:length(names)-1 names nnzr symm_spmv_perf_opt symm_spmv_perf_alpha],"|")
