import numpy as np

def append(f, data, name):
    
    ### get current shape 
    tmp = f[name].shape
    ### resize
    if len(tmp) == 1:
        f[name].resize((tmp[0] + data.shape[0],))
        f[name][tmp[0]:] = data
    elif len(tmp) == 2:
        f[name].resize((tmp[0], tmp[1] + 1))
        if len(data.shape) < 2:
            f[name][:, tmp[1]:] = data[:, np.newaxis]
        else:
            f[name][:, tmp[1]:] = data
    else:
        print("cannot append data to HDF5 with more than 2 dimensions", file=sys.stderr) 
        sys.exit(-1)


