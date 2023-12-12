import pandas as pd
import numpy as np
import os
os.system('cls')

file_path = 'C:/Users/migue/OneDrive/Ambiente de Trabalho/EEG stuff/IEEE project/EEG-Situ_VS_EEG-fMRI/stims.xlsx'
imgs_file = pd.read_excel(file_path)

# Perform random shuffles on the data
permutation = np.random.permutation(imgs_file.index)
shuffled_df = imgs_file.loc[permutation]

'''
This portion of code evaluates if there is space at the end of the stimuli array
to fit all the target stimuli with the specified TARGET_SEPERATION distance
In case there is not enough space, the target stimuli at the end of the array
will be moved directly to the begining of the array where the second part
of the algorithm will distribute them according to the valid distance between
target stimuli.
'''

FIRST_N_NONTARGET_TRIALS = 5
TARGET_SEPERATION = 3

invalid = True

df_size = len(shuffled_df) - 1
shuffled_df = shuffled_df.to_numpy()

# Loop will run until invalid target positions are removed
while invalid:
    invalid = False
    target_count = 0
    for it, (_, label) in enumerate(shuffled_df[::-1, :]):
        if label == "target":
            target_count += 1
            if (it + 1 - target_count) / TARGET_SEPERATION < 2.0:
                # Second target element (the invalid one) on the reversed array
                row_to_move = shuffled_df[df_size - it].copy()

                # Remove the row from its current position
                shuffled_df = np.delete(shuffled_df, df_size - it, axis=0)

                # Insert the row at the beginning
                shuffled_df = np.vstack([row_to_move, shuffled_df])

                invalid = True
                break
    print("Looping...")

tc = 0
for trial in shuffled_df:
    print(trial)
    if(trial[1] == "target"):
        tc += 1

print("\n======== After looping =======")
print("\nSize: ", len(shuffled_df))
print("Nº targets: ", tc)
print("Unique: ", len(np.unique(shuffled_df)))

'''
This second part is where the consecutive target stimuli will be redistributed according
to the valid target seperation distance. Furthermore, no target stimuli must be presented
in the first second of the trial
'''

non_target_count = 0
allocated_count = 0

allocated = []
reorganized_df = []

for it, (image, label) in enumerate(shuffled_df):
    if (it <= FIRST_N_NONTARGET_TRIALS and label == "non_target"): 
        reorganized_df.append([image, label])
        non_target_count += 1

    elif (it <= FIRST_N_NONTARGET_TRIALS and label == "target"):
        allocated.append([image, label])
        allocated_count += 1

    elif(it >= FIRST_N_NONTARGET_TRIALS):
        
        if (label == "target" and non_target_count >= TARGET_SEPERATION):
            reorganized_df.append([image, label])
            non_target_count = 0

        elif(label == "target" and non_target_count < TARGET_SEPERATION):
            
            allocated.append([image, label])
            allocated_count += 1

        elif(label == "non_target" and non_target_count >= TARGET_SEPERATION):
            if(allocated_count > 0):  
                reorganized_df.append(allocated[0])
                allocated_count -= 1
                allocated.pop(0) # To avoid duplicates
                
                non_target_count = 0
                reorganized_df.append([image, label])
                non_target_count += 1

            else:
                reorganized_df.append([image, label])
                non_target_count += 1

        elif(label == "non_target" and non_target_count < TARGET_SEPERATION):
            reorganized_df.append([image, label])
            non_target_count += 1

tc = 0
for trial in reorganized_df:
    print(trial[0])
    if(trial[1] == "target"):
        tc += 1

print("\n======== After distributing targets =======")
print("\nSize: ", len(reorganized_df))
print("Nº targets: ", tc)
print("Unique: ", len(np.unique(reorganized_df)))

print("\n======")
lol = [row[0] for row in reorganized_df]
print(lol)

reorganized_df = pd.DataFrame(reorganized_df, columns = ["image", "stim"])

# Check for duplicates based on entire rows
duplicates = reorganized_df[reorganized_df.duplicated(subset=["image"])]

if not duplicates.empty:
    print("Duplicates found based on entire rows:")
    print(duplicates)
else:
    print("No duplicates found.")

