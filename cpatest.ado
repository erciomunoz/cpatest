******************************
* cpatest.ado
* Ercio Munoz
* Graduate Center, CUNY
*! version 1.0  24-July-2017
*
* Computes Giacomini and White (2006 Econometrica) test
* It's based on Giacomini's matlab code
******************************

cap prog drop cpatest
program cpatest
version 14.2
syntax varlist(min=2 max=2 numeric) [if] [in] [, tau(integer 1) type(integer 1)]
marksample touse
local tau = `tau'
local type = `type'

if `type'<1 | `type'>2 {
	di as err "Type must be 1 (unconditional) or 2 (conditional)"
	exit
}

if `tau'<1 {
	di as err "The forecasting horizon must be a positive integer"
	exit
}

gettoken y1 y2 : varlist
tempvar lossdiff1 av_diff_loss
quiet: gen `lossdiff1' = `y1' - `y2'
quiet: mean `lossdiff1'
quiet: matrix `av_diff_loss' = e(b) 
if `av_diff_loss'[1,1]<0 {
local sign = "(-)"
}
else {
local sign = "(+)"
}

preserve
quiet: drop if `lossdiff1'==.
mata: compute("`lossdiff1'" , "`touse'" , `tau' , `type')
restore

display "(1) `y1' - `y2' = 0"
display ""
if `type'==1 {   
display as text "Choice of test   : Unconditional"
}
else {
display as text "Choice of test   : Conditional"
}
display as text "Forecast Horizon : " `tau'
display as text "Test-statistic   : " r(stat) "`sign'"
display as text "       P-value   : " r(pval)
end

mata
void function compute(string scalar varname, string scalar touse, real scalar tau, real scalar type) {
lossdiff1 = st_data(.,varname, touse)
TT = rows(lossdiff1)

if (type==1) {
instruments = J(TT,1,1)
lossdiff = lossdiff1
T = TT
}
else {
instruments = (J(TT-tau,1,1) , lossdiff1[1..TT-tau,1])
lossdiff = lossdiff1[tau+1..TT,1...]
T = TT-tau
}

/* Create the regressor matrix given by lossdiff*ht', where ht is the matrix of instruments */
cols_inst = cols(instruments)
rows_inst = rows(instruments)
reg = -999*J(rows_inst,cols_inst,1)
for (jj=1 ; jj<=cols_inst ; jj++) {
	reg[1...,jj] = instruments[1...,jj]:*lossdiff
} 

if (tau==1) {
/* calculate the test stat as nR^2 from the regression of one on lossdiff*ht */
ones2 = J(T,1,1)
resbeta = qrsolve(reg,ones2)
err = ones2-reg*resbeta
err2 = err:^2
r2 =  1-mean(err2)
teststat = T*r2
q = cols(reg)
pval = 1-chi2(q,abs(teststat))
st_numscalar("r(stat)",teststat)
st_numscalar("r(pval)",pval)
}
else {
regbar = mean(reg, 1)'
nlags = tau-1
q = cols(reg)
n = rows(reg)
reg = reg - J(n,1,1)*mean(reg)
gamma = -999*J(nlags,q,1)
samplevar = reg'*reg*n^-1
omegahat = samplevar
if (nlags>0) {
	for (i=1 ; i<=nlags ; i++) {
		reglag = (J(i,q,0) \ reg[1..n-i,1...])
		gamma = reg'*reglag+reglag'*reg
		gamma = gamma*n^-1
		nlags1 = nlags+1
		weights = 1-i*nlags1^-1
		omegahat = omegahat+weights*gamma 
	}
}
teststat = T*regbar'*invsym(omegahat)*regbar
pval = 1-chi2(q,abs(teststat))
st_numscalar("r(stat)",teststat)
st_numscalar("r(pval)",pval)
}
}
end
