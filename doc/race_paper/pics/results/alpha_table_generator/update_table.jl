skx=readdlm("../skx/data_symm_spmv/RLM/alphas.txt",',')
skx_name=skx[2:end,2]
skx_alphas=skx[2:end,4]
skx_nnz=skx[2:end,3]

ivy=readdlm("../ivy/data_symm_spmv/RLM/alphas.txt",',')
ivy_name=ivy[2:end,2]
ivy_alphas=ivy[2:end,4]


reqd_field=[3,5,7,9,11]
table=Any[]

for i in 1:length(skx_name)
	curr_mtx_name=skx_name[i]
	if(ivy_name[i] == curr_mtx_name)
		println(curr_mtx_name)
		row = Any[]
		push!(row, string("{",i,"}") )
		push!(row, string("& {", curr_mtx_name, "}") )
		push!(row, string("& ", 1/skx_nnz[i]) )
		push!(row, string("& ", skx_alphas[i]) )
		push!(row, string("& ", ivy_alphas[i]) )
		push!(row,string("\\\\"))
		push!(table, row)
	else
		println("Name mismatch")
	end
end

println(table)
writedlm("table.txt",  table)
				

