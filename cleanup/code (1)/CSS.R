# paper parameterization
x = seq(0,10,0.01)
y_min = 50
y_max = 100
lambda = 2
ic50 = 0.5

y1 = y_min + (y_max-y_min)/(1+(ic50/x)^lambda) # the paper definition
y2 = y_min + (y_max-y_min)/(1+exp(-lambda*(log(x)-log(ic50)))) # b = -lambda according to drm definition
plot(x,y1)
lines(x,y2)

# auc
c2=log10(max(x))
c1=log10(min(x[-1]))
m=log10(ic50)
auc1 = y_min*(c2-c1)+(y_max-y_min)*(1/lambda)*log10((1+10^(lambda*(c2-m)))/(1+10^(lambda*(c1-m)))) #auc according to the paper parameterization
auc1
auc1/(100*(c2-c1))
plot(x,y1,log="x")


a = y_max # corresponding to the d parameter (coef[3]) in drm
d = y_min # corresponding to the c parameter (coef[2]) in drm
c = log10(ic50) # corresponding to the log10(e) parameter in drm (log10(coef[4]))
b = lambda # if drm gives b (coef[1]), then b must be inversed to get the lambda, i.e lambda = -b
x2 = c2
x1 = c1
auc2 = (((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) # auc according to the drm parameterization
auc2

auc3 = (((((a-d)*log10(1+10^(b*(c-x2))))/(b*log10(10)))+a*x2)-((((a-d)*log10(1+10^(b*(c-x1))))/(b*log10(10)))+a*x1)) # base does not matter
auc3

# example
res = drm(y1~x,fct=LL.4())
coef = res$coefficients
coef
y3 = coef[2] + (coef[3]-coef[2])/(1+exp(coef[1]*(log(x)-log(coef[4])))) # b = -lambda according to drm definition
lines(x,y3,col="red")
y4 = coef[3] + (coef[2]-coef[3])/(1+exp(-coef[1]*(log(x)-log(coef[4])))) # b = -lambda according to drm definition
lines(x,y4,col="blue")

# the correct parameterization
a = coef[3]
d = coef[2]
c = log10(coef[4])
b = -coef[1]
x2 = c2
x1 = c1

css1 = (((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1))
css2 = (((((a-d)*log10(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1))

css # auc1=auc2=auc3=css1=css2

# another example
library(drc)
library(Brobdingnag)
tmp = read.table(file="clipboard",header=T)
res = drm(tmp$Perc_Inh~as.numeric(rownames(tmp)),fct=LL.4())
plot(res)
coef = res$coefficients
coef

a = coef[3]
d = coef[2]
c = log10(coef[4])
b = -coef[1]
x2 = log10(max(as.numeric(rownames(tmp))))
x1 = log10(min(as.numeric(rownames(tmp)[-1])))

css = (((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1))

css
