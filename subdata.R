# Function that produce a subset of data
# @param:
# data - full data set
# nsize - sample size in the target subset
# psize - feature size in the target subset
# causal - index of causal features
# seed - set seeds
# return:
# subset of data with specified dimension
# index of causal in the subset

subdata <- function(data,label,nsize, psize, causal, seed=1998){
  n = nrow(data); p = ncol(data)
  causal_prop = length(causal)/p
  set.seed(seed + 2)
  num = floor(causal_prop*psize)
  causal_sub = sample(causal, size = num)
  set.seed(seed + 4)
  sub_non_causal = sample((1:p)[-causal], size = (psize - num))
  sub_p = sort(c(causal_sub, sub_non_causal))
  set.seed(seed + 8)
  sub_n = sample(1:n, size = nsize)
  data_sub = data[sub_n,sub_p]
  label_sub = label[sub_n]
  return(list(data_sub = data_sub,
              label_sub = label_sub,
              causal_sub = sort(causal_sub)))
}
