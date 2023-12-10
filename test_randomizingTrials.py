import pandas as pd
import numpy as np

file_path = 'C:/Users/migue/OneDrive/Ambiente de Trabalho/EEG stuff/IEEE project/EEG-Situ_VS_EEG-fMRI/stims.xlsx'
imgs_file = pd.read_excel(file_path)

# Perform random shuffles on the data
permutation = np.random.permutation(imgs_file.index)
shuffled_df = imgs_file.loc[permutation]

FIRST_N_NONTARGET_TRIALS = 5
TARGET_SEPERATION = 3

invalid_count = 1

df_size = len(shuffled_df) - 1

'''
This portion of code evaluates if there is space at the end of the stimuli array
to fit all the target stimuli with the specified TARGET_SEPERATION distance
In case there is not enough space, the target stimuli at the end of the array
will be moved directly to the begining of the array where the second part
of the algorithm will distribute them according to the valid distance between
target stimuli.
'''

while (invalid_count > 0): # Loop will run until invalid target positions are removed
    invalid_count = 0
    target_count = 0
    for it, (image, label) in enumerate(shuffled_df.iloc[::-1].values):
        if (label == "target"):
            target_count += 1
        if(it > 0 and 
        label == "target" and 
        (it + 1 - target_count) / TARGET_SEPERATION < 2.0):
            
            # Second target element (the invalid one) on the reversed DataFrame
            row_to_move = shuffled_df.iloc[df_size - it].copy()
            shuffled_df = shuffled_df.drop(df_size - it) # Remove the row from its current position
            shuffled_df = pd.concat([row_to_move.to_frame().T, shuffled_df], ignore_index=True)
            invalid_count += 1
    print("Looping...")

''' 
This second part is where the consecutive target stimuli will be redistributed according
to the valid target seperation distance. Furthermore, no target stimuli must be presented
in the first second of the trial
'''

non_target_count = 0
allocated_count = 0

allocated = []
reorganized_df = []

for it, (image, label) in enumerate(shuffled_df.values):
    if (it < FIRST_N_NONTARGET_TRIALS and label == "non_target"): 
        reorganized_df.append([image, label])
        non_target_count += 1
        print("lol1", label )

    else:
        
        if (label == "target" and non_target_count >= TARGET_SEPERATION and it >= FIRST_N_NONTARGET_TRIALS):
            reorganized_df.append([image, label])
            non_target_count = 0
            print("lol2", label )

        elif(label == "target" and non_target_count < TARGET_SEPERATION):
            allocated.append([image, label])
            allocated_count += 1
            print("lol3", label )

        elif(label == "non_target" and non_target_count >= TARGET_SEPERATION and it >= FIRST_N_NONTARGET_TRIALS):
            if(allocated_count > 0):  
                reorganized_df.append(allocated[allocated_count - 1])
                allocated_count -= 1
                print("lol5", label)
                non_target_count = 0
                reorganized_df.append([image, label])
                non_target_count += 1
            else:
                reorganized_df.append([image, label])
                non_target_count += 1
                print("lol7", label)
        elif(label == "non_target" and non_target_count < TARGET_SEPERATION):
            reorganized_df.append([image, label])
            non_target_count += 1
            print("lol8", label)
        

print("\n----------") 
for trial in range(len(shuffled_df.values)):
    print(shuffled_df.values[trial][1])

print("\n----------")
targ_count = 0            
for trial in range(len(reorganized_df)):
    print(reorganized_df[trial][1])
    if(reorganized_df[trial][1] == "target"):
        targ_count += 1

print("\nSize before: ", len(shuffled_df))
print("Nº targets   : ", len(shuffled_df[shuffled_df["stimulus"] == "target"]))
unique = np.unique(shuffled_df.values)
print("Repeated elements: ", len(unique))

print("\nSize after : ", len(reorganized_df))
print("Nº targets   : ", targ_count)
unique = np.unique(reorganized_df)
print("Repeated elements: ", len(unique))

#print(np.sort(shuffled_df.values) + " || " + np.sort(reorganized_df))
