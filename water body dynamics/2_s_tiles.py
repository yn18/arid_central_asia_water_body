def read_tiles(tasks):
    from osgeo import gdal,ogr,osr

    tasks_ds = ogr.Open(tasks, 0)
    tasks_lyr = tasks_ds.GetLayer(0)

    tiles = []

    for tile_feat in tasks_lyr:
        tiles.append(tile_feat.GetField('tag'))
        
    del tasks_ds

    return tiles


def retiles(tasks, tile, tif, pathout):
    from osgeo import gdal,gdalconst,ogr,osr
    import os

    ds_tasks = ogr.Open(tasks, 0)
    lyr_tasks = ds_tasks.GetLayer(0)

    spatialRef_tasks = lyr_tasks.GetSpatialRef()
    spatialRef_tasks_wkt = spatialRef_tasks.ExportToWkt()# 投影信息

    feature_tile = lyr_tasks.GetNextFeature()
    geometry_tile = feature_tile.geometry().Clone()

    ring = geometry_tile.GetGeometryRef(0)
    point_1 = ring.GetPoints()[0]# tile左上角坐标
    point_3 = ring.GetPoints()[2]# tile右下角坐标

    expression = "{} = '{}'".format('tag', tile)
    ds_clip = gdal.Warp('', tif, format='MEM', cutlineDSName=tasks, cutlineWhere=expression, cropToCutline=True, outputBounds=(point_1[0],point_3[1],point_3[0],point_1[1]), outputBoundsSRS=spatialRef_tasks_wkt, xRes=0.00025, yRes=0.00025)

    pathout_1 = pathout + '/' + tile[0:4] + '/' + tile[-4:] + '/' + tile

    if os.path.isdir(pathout_1):
        pass
    else:
        try:
            os.makedirs(pathout_1)
        except:
            pass

    out_data = pathout_1 + '/' + tile + '_s' + '.tif'

    driver = gdal.GetDriverByName('GTiff')
    driver.CreateCopy(out_data, ds_clip, strict=1, options=["TILED=YES", "COMPRESS=LZW"])

    del ds_clip

    return


def divide(datas, n):
    '''进程分割'''

    mpi_datas = {}
    step = len(datas)//n
    for i in range(n):
        if i < n-1:
            mpi_data = datas[i*step:(i+1)*step]
            mpi_datas[i] = mpi_data
        else:
            mpi_data = datas[i*step:]
            mpi_datas[i] = mpi_data

    j = 0
    while len(mpi_datas[n-1]) > step and j < n-1:
        mpi_datas[j].append(mpi_datas[n-1][-1])
        mpi_datas[n-1].remove(mpi_datas[n-1][-1])
        j = j + 1
    
    mpi_datas_out = []
    for mpi_data_out in mpi_datas.values():
        mpi_datas_out.append(mpi_data_out)
    return mpi_datas_out


def main():
    import mpi4py.MPI as MPI
    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()

    import random
    import os

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ta', '--tasks', type=str, help='tasks', required=True)# tasks
    parser.add_argument('-i', '--input', type=str, help='input', required=True)# 输入路径
    parser.add_argument('-o', '--output', type=str, help='output', required=True)# 输出路径
    args = parser.parse_args()

    if comm_rank == 0:
        datas = read_tiles(args.tasks)
        random.shuffle(datas)
        mpi_datas = divide(datas, comm_size)
    else:
        datas = None
        mpi_datas = None
        tiles = None

    mpi_datas_divide = comm.scatter(mpi_datas, root=0)

    if os.path.isdir(args.output):
        pass
    else:
        try:
            os.makedirs(args.output)
        except:
            pass

    for tile in mpi_datas_divide:
        retiles(args.tasks, tile, args.input, args.output)
            
    return

if __name__ == '__main__':
    main()