import enum
import vtk
import vtk.util.numpy_support as vtknpsuppot
import numpy as np
import time
import csv

def rawImageReader(fname, nrows, ncols, nsls, dtype='uint8'):
    fd = open(fname)
    print('reading image data:',fname)
    t = time.time()
    A = np.fromfile(fd, dtype=dtype, sep="", count=nrows*ncols*nsls)
    fd.close()
    A = A.reshape([nsls,nrows,ncols])
    print('---> reading completed, time consumed: ', time.time()-t,'s')
    return A

def rawImageWriter(A, fname):
    fnameS = fname + '_' + \
            str(A.shape[2]) + 'x' + \
            str(A.shape[1]) + 'x' + \
            str(A.shape[0]) + '_' + \
            str(A.dtype) + '.raw'
    f = open(fnameS, "wb")
    print('saving raw image into ', fnameS)
    t = time.time()
    _ = f.write(A)
    f.close()
    print('  ----> saving complete (', str(int(time.time()-t)), 's )')

def vtuRead_old(fnameVTU):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fnameVTU)
    reader.Update()  # Needed because of GetScalarRange
    coords = reader.GetOutput().GetPoints().GetData()
    tMap = reader.GetOutput().GetPointData().GetArray("thickness")
    Nx = reader.GetOutput().GetPointData().GetArray("normal_x")
    Ny = reader.GetOutput().GetPointData().GetArray("normal_y")
    Nz = reader.GetOutput().GetPointData().GetArray("normal_z")
    coords = vtknpsuppot.vtk_to_numpy(coords)
    Nx = vtknpsuppot.vtk_to_numpy(Nx)
    Ny = vtknpsuppot.vtk_to_numpy(Ny)
    Nz = vtknpsuppot.vtk_to_numpy(Nz)
    tMap = vtknpsuppot.vtk_to_numpy(tMap)
    return coords, Nx, Ny, Nz, tMap

def vtuRead(fnameVTU, fieldnames=None):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fnameVTU)
    reader.Update()  # Needed because of GetScalarRange
    
    # read the coordinates of points
    coords = vtknpsuppot.vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
    
    # check the field names
    if fieldnames is None:
        fieldnames = []
        
        with open(fnameVTU, 'r') as f:
                        
            A = 'aaaaaaaaaa'
            while A != '<PointData':
                line = f.readline()                
                if len(line)>=10:
                    A = line[:10]
            
            line = f.readline()
            A = line[:10]
            while A == '<DataArray':
                id0 = line.find('Name="') + 6
                id1 = line.find('" NumberOfComponents')
                
                fieldnames.append(line[id0:id1])
                
                line = f.readline()
                A = line[:10]

    # read field data
    fields = {}
    for i, fieldname in enumerate(fieldnames):
        fields[fieldname] = vtknpsuppot.vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(fieldname))

    return coords, fields



'''
convert vtu unstructured grid to csv file
In the csv file, the first 3 colomns are x, y, z; and the following colomns are field data
'''
def vtu2csv(fname, fnameS=None, fieldnames=None, decimals=4):
    if fnameS is None:
        fnameS = fname[:-4] + '.csv'

    # load the vtu data
    print('reading vtu ...')
    coords, fields = vtuRead(fname, fieldnames)
    points = coords
    col_header = ['x', 'y', 'z']
    for key in fields:
        points = np.concatenate((points, np.array([fields[key]]).T), axis=1)
        col_header.append(key)

    # save into csv
    print('saving csv ...')
    with open(fnameS, 'w', newline='') as f:
        # create the csv writer
        writer = csv.writer(f, delimiter = ',')
        writer.writerow(col_header)
        writer.writerows(np.round(points, decimals)) #decimals - number of digits after .
    print('csv field saved: ' + fnameS)