import numpy as np
from scipy.cluster import vq

def python_vq(all_data,code_book):
    import time
    t1 = time.time()
    codes1,dist1 = vq.vq(all_data,code_book)
    t2 = time.time()
    #print 'fast (double):', t2 - t1
    #print '  first codes:', codes1[:5]
    #print '  first dist:', dist1[:5]
    #print '  last codes:', codes1[-5:]
    #print '  last dist:', dist1[-5:]
    float_obs = all_data.astype(np.float32)
    float_code = code_book.astype(np.float32)
    t1 = time.time()
    codes1,dist1 = vq.vq(float_obs,float_code)
    t2 = time.time()
    #print 'fast (float):', t2 - t1
    #print '  first codes:', codes1[:5]
    #print '  first dist:', dist1[:5]
    #print '  last codes:', codes1[-5:]
    #print '  last dist:', dist1[-5:]

    return codes1,dist1

def read_data(name):
    f = open(name,'r')
    data = []
    for line in f.readlines():
        data.append(map(float,string.split(line)))
    f.close()
    return array(data)

def main():
    np.random.seed((1000,1000))
    Ncodes = 40
    Nfeatures = 16
    Nobs = 4000
    code_book = np.random.normal(0,1,(Ncodes,Nfeatures))
    features = np.random.normal(0,1,(Nobs,Nfeatures))
    codes,dist = python_vq(features,code_book)

if __name__ == '__main__':
    main()
