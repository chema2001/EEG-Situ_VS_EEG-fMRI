import os
import numpy as np

def run_fast_scandir(dir, ext):    # dir: str, ext: list
    subfolders, datafiles = [], []

    for f in os.scandir(dir):
        if f.is_dir():
            subfolders.append(f.path)
        if f.is_file():
            if os.path.splitext(f.name)[1].lower() in ext:
                datafiles.append(f.path)

    for dir in list(subfolders):
        sf, f = run_fast_scandir(dir, ext)
        subfolders.extend(sf)
        datafiles.extend(f)
    return subfolders, datafiles

_, imgs = run_fast_scandir('C:/Users/migue/OneDrive/Ambiente de Trabalho/EEG stuff/IEEE project/EEG-Situ_VS_EEG-fMRI/Celebs', [".jpg"])

target_img = 'Daniel_Radcliffe'

from pandas import DataFrame

images = []
stim = []

for i, img in enumerate(imgs):
    aux = img.split('\\')[:]
    conc_string = ''.join(["Celebs"+"/", aux[1]])
    images.append(conc_string)
    print(conc_string)

    if target_img in conc_string:
        aux_stim = 'target'
    else:
        aux_stim = 'non_target'
    stim.append(aux_stim)
    
df = DataFrame({'image': images, 'stimulus': stim})

'''
permutation = np.random.permutation(df.index)
shuffled_df = df.loc[permutation] # .loc -> Select rows
shuffled_df.to_excel('stims.xlsx', sheet_name='sheet1', index=False)
'''