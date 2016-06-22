

rootComb = function(x,root = 1){ 
  n  = sum(x)
  co = 1
  for(i in x){
    #cat(n,' ',i,'\n')
    co = co * rootChoose(n,i,root)
    n = n - i
    }
  return(co)
  }


stabilityExp = function(x,root = 1){
  co  = 0
  n  = sum(x)
  totComb = rootComb(x,root)
  for(i in 1:length(x)){
    k = x[i]
    if(k > 1){
      rootCo     = rootChoose(n-2,k-2,root) * rootComb(x[-i],root)
      normRootCo = rootCo / totComb
      normCo     = normRootCo ** root
      co = co + normCo
      } else {
      rootCo = rootChoose(n-2,k-2,root) * rootComb(x[-i],root)
      normRootCo = rootCo / totComb
      normCo = normRootCo ** root
      co = co + normCo     
      }
    }
  return(co)
  }


instab = function(x,y,root){
  a = stabilityExp(x,root)
  b = stabilityExp(y,root)
  return(a * (1-b) + (1-a) * b)
  }

# 
# stabilityExp(c(100,10,1,1),10)
# 
# rootChoose(10,2,1)*rootChoose(10,2,1)
# 
# a = stabilityExp(c(1,1,1,3,5))
# b = stabilityExp(c(1,1,1,3,5))
# a * (1-b) + (1-a) * b
# 
# 
# rootComb(c(1,1,2,2),1)
# 
# choose(6,1) * choose(5,1) * choose(4,2) * choose(2,2)








