from matplotlib.pyplot import savefig, subplots, show, title, suptitle, rc 
from numpy import array, zeros, arange
#rc('text', usetex = True)
nodes_number = 4

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

lines = 0  
section_witdh = nodes_number * 2
data = []

with open('data/times.dat') as file:	
	for line in file.readlines()[1:]:
		data += [[num(x) for x in line.split()[1:]]]
		lines += 1

sections = lines // section_witdh

for s in range(sections):
	homemade = data[s * section_witdh : nodes_number + s * section_witdh]
	fftw_mpi = data[nodes_number + s * section_witdh : 2 * nodes_number + s * section_witdh]
	versions = (homemade, fftw_mpi)

	nx = homemade[0][1]
	ny = homemade[0][2]
	nz = homemade[0][3]
	itrs = homemade[0][4]
	
	total = [array([t for *_, t in version]) for version in versions]
	
	nodes = array([2 ** i for i in range(nodes_number)])

	fig, ax = subplots(layout = 'constrained')
	bottom = zeros(nodes_number)
	width = 0.2
	i = 0
	multiplier = 0
	bar_colors = ('red', 'yellow')
	label = ('Homemade', 'Western')


	for times in total:
		offset = width * multiplier - width / 2

		p = ax.bar(nodes + offset, times, width, label = label[i], bottom = bottom, \
					     color = bar_colors[i], edgecolor = 'black')
					
		multiplier += 1
		i += 1     

	ax.set_xticks(nodes)
	ax.set_ylabel(r'Time [$s$]')
	ax.set_xlabel('Number of nodes')
	ax.legend(loc = 'upper right')
	ax.grid(linestyle = '--', axis = 'y')

	# rc('text', usetex=True)
	title(r'$\bf{Total}$ $\bf{time}$ $\bf{per}$ $\bf{number}$ $\bf{of}$ $\bf{nodes}$' + \
				f'\nSize of grid: {nx}' + r'$\times$' + f'{ny}' + r'$\times$' + f'{nz}' + \
				f'\nIterations: {itrs}', fontsize = 10)

	savefig(f'analysis/{nx}_{ny}_{nz}_{itrs}.png')
