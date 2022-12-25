import os
import rasterio
import numpy as np
from tqdm import tqdm
import time
import pandas as pd


def read_tiles(tasks):
    from osgeo import gdal,ogr,osr

    tasks_ds = ogr.Open(tasks, 0)
    tasks_lyr = tasks_ds.GetLayer(0)

    tiles = []

    for tile_feat in tasks_lyr:
        tiles.append(tile_feat.GetField('tag'))
        
    del tasks_ds

    return tiles


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


def code2idx(code, imwidth):
    row = int(code / imwidth)
    column = code - row * imwidth
    return [row, column]


def inarea(idx, imsize):
    if 0 <= idx[0] < imsize[0] and 0 <= idx[1] < imsize[1]:
        return True
    else:
        return False


def dfs(idx, imsize, lbl, hydro_band, jrc_band, save_band, water_code):
    stack = []
    hydro_band[idx[0], idx[1]] = 0
    stack.append(idx)

    stats = [0, 0]
    while stack:
        idx = stack[-1]
        if inarea([idx[0] - 1, idx[1]], hydro_band.shape) and hydro_band[idx[0] - 1, idx[1]] == lbl:
            idx = [idx[0] - 1, idx[1]]
            hydro_band[idx[0], idx[1]] = 0
            stack.append(idx)
        elif inarea([idx[0], idx[1] + 1], hydro_band.shape) and hydro_band[idx[0], idx[1] + 1] == lbl:
            idx = [idx[0], idx[1] + 1]
            hydro_band[idx[0], idx[1]] = 0
            stack.append(idx)
        elif inarea([idx[0] + 1, idx[1]], hydro_band.shape) and hydro_band[idx[0] + 1, idx[1]] == lbl:
            idx = [idx[0] + 1, idx[1]]
            hydro_band[idx[0], idx[1]] = 0
            stack.append(idx)
        elif inarea([idx[0], idx[1] - 1], hydro_band.shape) and hydro_band[idx[0], idx[1] - 1] == lbl:
            idx = [idx[0], idx[1] - 1]
            hydro_band[idx[0], idx[1]] = 0
            stack.append(idx)
        else:
            idx = stack[-1]
            if_nodata, stats_1 = mark(idx, imsize, lbl, jrc_band, save_band, water_code)
            stats[0] += stats_1[0]
            stats[1] += stats_1[1]
            # if if_nodata:
            #     save_band[save_band == lbl] = 0
            #     return
            stack.pop()
    return stats


def mark(idx, imsize, lbl, jrc_band, save_band, water_code):
    stated = 4
    nodata = 0
    stats = [0, 0]
    if jrc_band[idx[0], idx[1]] == nodata:
        jrc_band[idx[0], idx[1]] = stated
        stats[1] += 1
        # return True
    if inarea([idx[0] - 1, idx[1]], imsize) and jrc_band[idx[0] - 1, idx[1]] == nodata:
        jrc_band[idx[0] - 1, idx[1]] = stated
        stats[1] += 1
        # return True
    if inarea([idx[0] - 1, idx[1] + 1], imsize) and jrc_band[idx[0] - 1, idx[1] + 1] == nodata:
        jrc_band[idx[0] - 1, idx[1] + 1] = stated
        stats[1] += 1
        # return True
    if inarea([idx[0], idx[1] + 1], imsize) and jrc_band[idx[0], idx[1] + 1] == nodata:
        jrc_band[idx[0], idx[1] + 1] = stated
        stats[1] += 1
        # return True
    if inarea([idx[0] + 1, idx[1] + 1], imsize) and jrc_band[idx[0] + 1, idx[1] + 1] == nodata:
        jrc_band[idx[0] + 1, idx[1] + 1] = stated
        stats[1] += 1
        # return True
    if inarea([idx[0] + 1, idx[1]], imsize) and jrc_band[idx[0] + 1, idx[1]] == nodata:
        jrc_band[idx[0] + 1, idx[1]] = stated
        stats[1] += 1
        # return True
    if inarea([idx[0] + 1, idx[1] - 1], imsize) and jrc_band[idx[0] + 1, idx[1] - 1] == nodata:
        jrc_band[idx[0] + 1, idx[1] - 1] = stated
        stats[1] += 1
        # return True
    if inarea([idx[0], idx[1] - 1], imsize) and jrc_band[idx[0], idx[1] - 1] == nodata:
        jrc_band[idx[0], idx[1] - 1] = stated
        stats[1] += 1
        # return True
    if inarea([idx[0] - 1, idx[1] - 1], imsize) and jrc_band[idx[0] - 1, idx[1] - 1] == nodata:
        jrc_band[idx[0] - 1, idx[1] - 1] = stated
        stats[1] += 1
        # return True

    if jrc_band[idx[0], idx[1]] == water_code:
        jrc_band[idx[0], idx[1]] = stated
        stats[0] += 1
        save_band[idx[0], idx[1]] = lbl
    if inarea([idx[0] - 1, idx[1]], imsize) and jrc_band[idx[0] - 1, idx[1]] == water_code:
        jrc_band[idx[0] - 1, idx[1]] = stated
        stats[0] += 1
        save_band[idx[0] - 1, idx[1]] = lbl
    if inarea([idx[0] - 1, idx[1] + 1], imsize) and jrc_band[idx[0] - 1, idx[1] + 1] == water_code:
        jrc_band[idx[0] - 1, idx[1] + 1] = stated
        stats[0] += 1
        save_band[idx[0] - 1, idx[1] + 1] = lbl
    if inarea([idx[0], idx[1] + 1], imsize) and jrc_band[idx[0], idx[1] + 1] == water_code:
        jrc_band[idx[0], idx[1] + 1] = stated
        stats[0] += 1
        save_band[idx[0], idx[1] + 1] = lbl
    if inarea([idx[0] + 1, idx[1] + 1], imsize) and jrc_band[idx[0] + 1, idx[1] + 1] == water_code:
        jrc_band[idx[0] + 1, idx[1] + 1] = stated
        stats[0] += 1
        save_band[idx[0] + 1, idx[1] + 1] = lbl
    if inarea([idx[0] + 1, idx[1]], imsize) and jrc_band[idx[0] + 1, idx[1]] == water_code:
        jrc_band[idx[0] + 1, idx[1]] = stated
        stats[0] += 1
        save_band[idx[0] + 1, idx[1]] = lbl
    if inarea([idx[0] + 1, idx[1] - 1], imsize) and jrc_band[idx[0] + 1, idx[1] - 1] == water_code:
        jrc_band[idx[0] + 1, idx[1] - 1] = stated
        stats[0] += 1
        save_band[idx[0] + 1, idx[1] - 1] = lbl
    if inarea([idx[0], idx[1] - 1], imsize) and jrc_band[idx[0], idx[1] - 1] == water_code:
        jrc_band[idx[0], idx[1] - 1] = stated
        stats[0] += 1
        save_band[idx[0], idx[1] - 1] = lbl
    if inarea([idx[0] - 1, idx[1] - 1], imsize) and jrc_band[idx[0] - 1, idx[1] - 1] == water_code:
        jrc_band[idx[0] - 1, idx[1] - 1] = stated
        stats[0] += 1
        save_band[idx[0] - 1, idx[1] - 1] = lbl
    return False, stats


def extract_lake(hydro_path, jrc_path, save_path, year, water_code):
    hydro_ds = rasterio.open(hydro_path)
    jrc_ds = rasterio.open(jrc_path)

    hydro_band = hydro_ds.read(1)
    jrc_band = jrc_ds.read(1)

    hydro_size = hydro_band.shape
    unique_values = [x for x in np.unique(hydro_band) if x != hydro_ds.nodata]
    unique_counts = [np.sum(hydro_band == x) for x in unique_values]
    unique_seeds = [np.random.choice(np.where(hydro_band.flat == x)[0]) for x in unique_values]
    uniques = {a: [b, c] for a, b, c in zip(unique_values, unique_counts, unique_seeds)}
    uniques = sorted(uniques.items(), key=lambda x: x[1][0], reverse=True)

    save_band = np.zeros((hydro_size[0], hydro_size[1]), dtype=hydro_band.dtype)
    stats_all = {'label': [], 'lake': [], 'nodata': [], 'year': []}
    for unique in uniques:
        hydro_lbl = unique[0]
        seed = unique[1][1]
        idx = code2idx(seed, hydro_size[1])
        stats_1 = dfs(idx, hydro_band.shape, hydro_lbl, hydro_band, jrc_band, save_band, water_code)
        stats_all['label'].append(hydro_lbl)
        stats_all['lake'].append(stats_1[0])
        stats_all['nodata'].append(stats_1[1])
        stats_all['year'].append(year)

    dst_ds = rasterio.open(save_path, 'w', driver='GTiff', height=hydro_size[0], width=hydro_size[1], count=1,
                           dtype=hydro_band.dtype, crs=jrc_ds.crs, transform=jrc_ds.transform)
    dst_ds.nodata = 0
    dst_ds.write(save_band, 1)

    return stats_all


def main():
    import mpi4py.MPI as MPI
    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()

    from osgeo import gdal,ogr,osr
    import pandas as pd
    import os
    import random
    import shutil

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ihy', '--input_h', type=str, help='input_h', required=True)# hydro_path
    parser.add_argument('-ijr', '--input_j', type=str, help='input_j', required=True)# jrc_path
    parser.add_argument('-ta', '--tasks', type=str, help='tasks', required=True)# task_shp
    parser.add_argument('-o', '--output', type=str, help='output', required=True)# mask
    args = parser.parse_args()

    if comm_rank == 0:
        datas = read_tiles(args.tasks)
        random.shuffle(datas)
        mpi_datas = divide(datas, comm_size)
        years = [x+1998 for x in range(23)]

        if os.path.exists(args.output):
            shutil.rmtree(args.output)
        os.mkdir(args.output)
    else:
        datas = None
        mpi_datas = None
        years = None

    mpi_data_divide = comm.scatter(mpi_datas, root=0)
    years = comm.bcast(years, root=0)

    stats_all = {'label': [], 'lake': [], 'nodata': [], 'year': []}
    for tile in mpi_data_divide:
        for year in years:

            hydro_path = args.input_h + '/' + tile[0:4] + '/' + tile[-4:] + '/' + tile + '/' + tile + '_s.tif'
            jrc_path = args.input_j + '/' + tile[0:4] + '/' + tile[-4:] + '/' + tile + '/' + tile + '_' + str(year) + '.tif'
            save_dir = args.output + '/' + tile[0:4] + '/' + tile[-4:] + '/' + tile
            save_path = save_dir + '/' + tile + '_' + str(year) + '_o.tif'

            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            stats_1 = extract_lake(hydro_path, jrc_path, save_path, year, 3)
            stats_all['label'].extend(stats_1['label'])
            stats_all['lake'].extend(stats_1['lake'])
            stats_all['nodata'].extend(stats_1['nodata'])
            stats_all['year'].extend(stats_1['year'])
        print(tile)

    with open('/public/home/mfeng/ynsu/1_1/data/stats.txt', 'w') as f:
        f.write('label, lake, nodata, year\n')
        for i in range(len(stats_all['label'])):
            f.write('{}, {}, {}, {}\n'.format(stats_all['label'][i], stats_all['lake'][i], stats_all['nodata'][i], stats_all['year'][i]))
    # stats_path = args.output + '/' + 'stats.xlsx'
    # dst_df = pd.DataFrame(stats_all)
    # dst_df.to_excel(stats_path)

    return


if __name__ == "__main__":
    main()

