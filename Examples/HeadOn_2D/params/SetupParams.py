import os
import copy
import numpy as np

output_dir = './v1_withGW/'

q_array = np.linspace(1,10,19)

#Fa = 0.5

config = open(output_dir+'template.txt','r')
config_lines = []
for line in config.readlines():
  config_lines.append(line)

for q in q_array:
    
    temp_config = copy.deepcopy(config_lines)
    index_str = str(q).replace('.','d')
#    try:
#        os.mkdir('/home/dc-wong1/data/dc-wong1/BBH/HeadOn/2D/v1_withGW/q'+index_str)
#    except:
#        continue
    
    temp_config[1] = 'chk_prefix = /home/dc-wong1/data/dc-wong1/BBH/HeadOn/2D/v1_withGW/q'+index_str+'/HeadOn2D_\n'
    temp_config[2] = 'plot_prefix = /home/dc-wong1/data/dc-wong1/BBH/HeadOn/2D/v1_withGW/q'+index_str+'/HeadOn2DPlot_\n'
    temp_config[3] = 'extraction_prefix = /home/dc-wong1/data/dc-wong1/BBH/HeadOn/2D/v1_withGW/q'+index_str+'/HeadOn2D_\n'
    temp_config[4] = 'restart_file = /home/dc-wong1/data/dc-wong1/BBH/HeadOn/2D/v1_withGW/q'+index_str+'/HeadOn2D_002000.2d.hdf5\n'
    temp_config[130] = 'massB = '+str(q*0.5)+'\n'
    temp_config[128] = 'offsetA = '+str(q*10.0)+' 0.0\n'
    temp_output = open(output_dir+'/q'+index_str+'.txt','w')
    temp_output.writelines(temp_config)
    temp_output.close()

