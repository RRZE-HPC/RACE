using LsqFit

mkl=[
  40.0   64000.0          9.92   968724.0      
  60.0  216000.0         76.31        2.73946e6
  70.0  343000.0        230.26        5.06638e6
  80.0  512000.0        452.685       7.67008e6
  90.0  729000.0       1068.71        1.21006e7
 120.0       1.728e6   6600.29        3.61523e7
 140.0       2.744e6  16473.6         6.5312e7 
]

cgmn=[
40		64000	78.99	332364
60		216000	1007  661684
70		343000	 2277.09  959544
80		512000 4302.18 1440628
90		729000 7839 1827276
120	1728000 29777.8 4106652
140 	2744000 57916.9 6316176
#160	4096000 118008 9379604
#220	10648000 527629.5 23632788
]


feast_model_mkl(x, p) = 76.31*(x.^3./(60^3)).^p[1]
feast_model_cgmn(x, p) = 1007*(x.^3./(60^3)).^p[1]

feast_model_mkl_mem(x, p) = 2.73946e6*(x.^3./(60^3)).^p[1]
feast_model_cgmn_mem(x, p) = 661684*(x.^3./(60^3)).^p[1]
#feast_model(x, p) = p[1]*exp.(-x.*p[2])

function fit_param_mkl()
	xdata=mkl[2:end,1]
	ydata=mkl[2:end,3]
	p0=[0.5]
	println(xdata)
	println(ydata)
	curve_fit(feast_model_mkl, xdata, ydata, p0)
end

function fit_param_cgmn()
	xdata=cgmn[2:end,1]
	ydata=cgmn[2:end,3]
	p0=[0.5]
	println(xdata)
	println(ydata)
	curve_fit(feast_model_cgmn, xdata, ydata, p0)
end

function fit_param_mkl_mem()
	xdata=mkl[2:end,1]
	ydata=mkl[2:end,4]
	p0=[0.5]
	println(xdata)
	println(ydata)
	curve_fit(feast_model_mkl_mem, xdata, ydata, p0)
end

function fit_param_cgmn_mem()
	xdata=cgmn[2:end,1]
	ydata=cgmn[2:end,4]
	p0=[0.5]
	println(xdata)
	println(ydata)
	curve_fit(feast_model_cgmn_mem, xdata, ydata, p0)
end

