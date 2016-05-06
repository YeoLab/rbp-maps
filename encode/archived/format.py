from scipy.stats import norm

import pandas as pd
x = pd.read_table('test2/excluded.txt')
print norm.pdf(x)
