import pandas as pd
import numpy as np
df = pd.read_csv('NO_diff.csv')
column = df['base_num_fraction']

# 检查该列中是否有空值
has_nan = column.isnull().any()
print(has_nan) 
#print(np.min(wga_df["g_value"]))
#print(len(wga_df['position'][wga_df["g_value"]<=0.1]))