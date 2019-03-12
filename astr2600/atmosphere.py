import pandas as pd
import pkg_resources
path = pkg_resources.resource_filename(__name__, 'data/earth-standard-atmosphere.txt')
table = pd.read_csv(path, sep='\s+')

altitude = table['altitude'].values
temperature = table['temperature'].values
