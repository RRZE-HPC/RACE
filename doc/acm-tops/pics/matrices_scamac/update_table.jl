mtx=readdlm("matrix.txt")
mtx_name=mtx[2:end,3]

verbose=mtx[2:end,:]

#strip .mtx at end of mtx_name
for i in 1:size(verbose)[1]
	start_idx=-1
	if(searchindex(verbose[i,3],".mtx")!=0)
		start_idx = search(verbose[i,3],".mtx")[1]
	else
		start_idx = search(verbose[i,3],".mm")[1]
	end
	name = verbose[i,3][1:start_idx-1]
	verbose[i,3] = name
end


reqd_field=[3,5,7,9,11]
table=Any[]

for i in 1:length(mtx_name)
	curr_mtx_name=mtx_name[i]
	for j in 1:length(verbose[:,1])
		if((string(verbose[j,3],".mtx") == curr_mtx_name) || (string(verbose[j,3],".mm") == curr_mtx_name) )
			println(curr_mtx_name)
			row = Any[]
			push!(row, string("{",i,"}"))
			k=reqd_field[1]
			push!(row, string("& {",verbose[j,k],"}"))
			for l in 2:length(reqd_field)
				k=reqd_field[l]
				push!(row, string("& ",verbose[j,k],""))
			end
			push!(row, string("& {} &"))
			push!(row,string("\\\\"))
			push!(table, row)
			break
		end
	end

end

println(table)
writedlm("table.txt",  table)
				

