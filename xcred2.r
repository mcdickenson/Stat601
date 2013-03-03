# Stat 601
# Matt Dickenson
# Extra credit problem

require('zipfR') # for Ibeta
require('stats') # for optim

betacize <- function(shapes, c=C, d=D, e=E, f=F){

	a = shapes[1]
	b = shapes[2]

	c.err = .05 - Ibeta(c, a, b)
	d.err = .95 - Ibeta(d, a, b)
	e.err = .025 - Ibeta(e, a, b)
	f.err = .975 - Ibeta(f, a, b)

	sq.err = c(c.err^2, d.err^2, e.err^2, f.err^2)
	ssq.err = sum(sq.err)
	return(ssq.err)

	# maxsq.err = max(sq.err)
	# return(maxsq.err)

	# abs.err = c(abs(c.err), abs(d.err), abs(e.err), abs(f.err))
	# sabs.err = sum(abs.err)
	# return(sabs.err)
}

start = c(1,1) # sets initial dist to uniform


# Example 1
C = .05
D = .95
E = .025
F = .975
out = optim(par=start, fn=betacize, method="SANN")
print(out$par)

qbeta(.05, 1,1)

# Example 2
C = .5
D = .9
E = .4
F = .95

out2 = optim(par=start, fn=betacize)
print(out2$par)

qbeta(.05, out2$par[1], out2$par[2])
qbeta(.95, out2$par[1], out2$par[2])
qbeta(.025, out2$par[1], out2$par[2])
qbeta(.975, out2$par[1], out2$par[2])


# Correctness test 
C = pbeta(.05, 1, 2)
D = pbeta(.95, 1, 2)
E = pbeta(.025, 1, 2)
F = pbeta(.975, 1, 2)
out = optim(par=start, fn=betacize)
print(out$par)
