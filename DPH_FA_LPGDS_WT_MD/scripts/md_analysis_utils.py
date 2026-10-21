"""Read OpenMM DCD coordinates and fit protein coordinates without side effects."""
import struct
import numpy as np
def read_dcd(path):
    def record(f):
        size=f.read(4)
        if not size:return None
        n=struct.unpack('<i',size)[0];data=f.read(n)
        if len(data)!=n or f.read(4)!=size:raise ValueError('Truncated DCD record')
        return data
    with open(path,'rb') as f:
        h=record(f)
        if h[:4]!=b'CORD':raise ValueError('Not a coordinate DCD')
        expected=struct.unpack('<i',h[4:8])[0];boxflag=struct.unpack('<i',h[44:48])[0];record(f);natoms=struct.unpack('<i',record(f))[0]
        xyz=[];boxes=[]
        for _ in range(expected):
            if boxflag!=1:raise ValueError('Missing unit cell')
            cell=np.frombuffer(record(f),dtype='<f8')
            if not np.allclose(cell[[1,3,4]],0,atol=1e-6):raise ValueError('Triclinic cell unsupported')
            x=np.column_stack([np.frombuffer(record(f),dtype='<f4') for i in range(3)])
            if x.shape!=(natoms,3):raise ValueError('DCD atom count mismatch')
            xyz.append(x);boxes.append(cell[[0,2,5]])
        if record(f) is not None:raise ValueError('DCD contains extra data')
    return np.asarray(xyz),np.asarray(boxes)
def kabsch(moving,target):
    mc=moving.mean(axis=0);tc=target.mean(axis=0)
    u,s,vt=np.linalg.svd((moving-mc).T@(target-tc))
    rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    return rotation,mc,tc
