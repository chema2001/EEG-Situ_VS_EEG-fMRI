import pandas as pd
import random
import numpy as np

file_path = 'C:/Users/migue/OneDrive/Ambiente de Trabalho/IEEE project/EEG-Situ_VS_EEG-fMRI/stims.xlsx'
imgs_file = pd.read_excel(file_path)

permutation = np.random.permutation(imgs_file.index)
shuffled_df = imgs_file.loc[permutation]

x = shuffled_df['image'].values

for trial in range(len(x)):
    print(x[trial])


from PIL import Image
im = Image.open(x[0])
im.show()
