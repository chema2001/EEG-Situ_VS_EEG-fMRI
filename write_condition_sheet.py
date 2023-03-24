import os

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


subfolders, datafiles = run_fast_scandir('C:/Users/migue/OneDrive/Ambiente de Trabalho/IEEE project/EEG-Situ_VS_EEG-fMRI/Dataset', [".jpg"])

from pandas import DataFrame

images=[]

for image in datafiles:
    aux = image.split('\\')[:]
    conc_string = ''.join(["Dataset/", aux[:][1]])
    print(conc_string)
    images.append(conc_string)
    
df = DataFrame({'stimulus_time': images})
df.to_excel('stims.xlsx', sheet_name='sheet1', index=False)
