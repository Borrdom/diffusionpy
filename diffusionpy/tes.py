
ncomp=4
def mat(vec):
    k=0
    m=np.zeros((ncomp,ncomp))
    for i in range(ncomp):
        for j in range(i):
            m[i, j] = vec[k]
            m[j, i] = vec[k]
            k+=1
    return m
vec=[1,2,3,4,5,6]
print(mat(vec))





ncomp=3
def mat(vec):
    m=np.zeros((ncomp,ncomp))
    for i,j,k in zip(*np.triu_indices(ncomp,k=1),range(len(vec)+1)): m[i,j]=vec[k];m[j,i]=vec[k]
    return m
vec=[1,2,3]
print(mat(vec))
